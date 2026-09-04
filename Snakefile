configfile: "config.yaml"

DESIGNS = []
EPITOPE_OF = {}
for ep_name, ep in config["epitopes"].items():
    for i in range(ep["n_designs"]):
        d = f"{ep_name}_design_{i}"
        DESIGNS.append(d)
        EPITOPE_OF[d] = ep_name

VARIANTS = list(range(1, config["rosetta"]["nstruct"] + 1))

rule all:
    input:
        "results/final_results_table.csv"

rule rfdiffusion:
    output:
        pdb = "results/00_rfdiffusion/{design}.pdb"
    params:
        env = config["envs"]["rfdiffusion"],
        target = config["target_pdb"],
        contigmap = config["contigmap"],
        hotspots = lambda wc: config["epitopes"][EPITOPE_OF[wc.design]]["hotspot_res"],
        prefix = lambda wc: f"results/00_rfdiffusion/{wc.design}",
        hydra_dir = lambda wc: f"results/00_rfdiffusion/hydra_logs/{wc.design}"
    shell:
        """
        mkdir -p {params.hydra_dir}
        GPU_ID=$(( RANDOM % 8 ))
        CUDA_VISIBLE_DEVICES=$GPU_ID conda run -n {params.env} --no-capture-output \
        python /home/apdivek/RFdiffusion/scripts/run_inference.py \
            hydra.run.dir={params.hydra_dir} \
            inference.input_pdb={params.target} \
            inference.output_prefix={params.prefix} \
            'contigmap.contigs={params.contigmap}' \
            'ppi.hotspot_res=[{params.hotspots}]' \
            inference.num_designs=1 \
            denoiser.noise_scale_ca=0 \
            denoiser.noise_scale_frame=0
        mv {params.prefix}_0.pdb {output.pdb}
        mv {params.prefix}_0.trb {output.pdb}.trb 2>/dev/null || true
        """

rule make_resfile:
    input:
        pdb = "results/00_rfdiffusion/{design}.pdb"
    output:
        resfile = "results/01_rosetta/{design}.resfile"
    params:
        env = config["envs"]["biopython"]
    shell:
        """
        conda run -n {params.env} --no-capture-output \
        python scripts/make_resfile.py {input.pdb} {output.resfile}
        """

rule rosetta_fixbb:
    input:
        pdb = "results/00_rfdiffusion/{design}.pdb",
        resfile = "results/01_rosetta/{design}.resfile"
    output:
        clean = "results/01_rosetta/{design}_000{variant}_clean.pdb"
    params:
        binary = config["rosetta"]["binary"],
        ex_flags = config["rosetta"]["ex_flags"],
        outdir = "results/01_rosetta/",
        variant_prefix = lambda wc: f"{wc.design}_v{wc.variant}_"
    shell:
        """
        {params.binary} \
            -in:file:s {input.pdb} \
            -resfile {input.resfile} \
            -out:prefix {params.variant_prefix} \
            -out:path:all {params.outdir} \
            -nstruct 1 \
            {params.ex_flags} \
            -overwrite
        RAW=$(find {params.outdir} -name "{params.variant_prefix}*0001.pdb" | head -1)
        if [ -z "$RAW" ]; then
            echo "ERROR: no Rosetta output found matching {params.variant_prefix}*0001.pdb"
            exit 1
        fi
        cp "$RAW" {output.clean}
        """

rule renumber:
    input:
        clean = "results/01_rosetta/{design}_000{variant}_clean.pdb"
    output:
        renumbered = "results/02_renumbered/{design}_000{variant}.pdb"
    params:
        env = config["envs"]["biopython"]
    shell:
        """
        conda run -n {params.env} --no-capture-output \
        python scripts/renumber.py {input.clean} {output.renumbered}
        """

rule af2_initial_guess:
    input:
        pdb = "results/02_renumbered/{design}_000{variant}.pdb"
    output:
        pred = "results/03_af2/{design}_000{variant}_af2pred.pdb",
        score = "results/03_af2/{design}_000{variant}_scores.sc"
    params:
        env = config["envs"]["af2"],
        script = config["af2"]["script"],
        recycle = config["af2"]["recycle"],
        indir = "results/02_renumbered/",
        outdir = "results/03_af2/",
        runlist_file = "results/03_af2/{design}_000{variant}.runlist"
    shell:
        """
        echo "{wildcards.design}_000{wildcards.variant}" > {params.runlist_file}
        conda run -n {params.env} --no-capture-output \
        python {params.script} \
            -pdbdir {params.indir} \
            -outpdbdir {params.outdir} \
            -scorefilename {output.score} \
            -runlist {params.runlist_file} \
            -recycle {params.recycle}
        """

rule aggregate:
    input:
        expand("results/03_af2/{design}_000{variant}_af2pred.pdb",
               design=DESIGNS, variant=VARIANTS)
    output:
        "results/final_results_table.csv"
    shell:
        "python scripts/aggregate_results.py results/03_af2/ {output}"
