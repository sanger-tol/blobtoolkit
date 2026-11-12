process BUSCO_PREPARE_PRECOMPUTED {
    tag "${meta.id}"
    label 'process_single'
    // Technically only Python is needed but using the blobtoolkit container since it's already pulled by the pipeline
    container "docker.io/genomehubs/blobtoolkit:4.4.6"

    input:
    tuple val(meta), path(multi_copy), path(single_copy), path(fragmented)

    output:
    tuple val(meta), path("busco_sequences.tar.gz"), emit: busco

    when: task.ext.when == null || task.ext.when

    script:
    """
    python - <<'PY'
    import tarfile
    from pathlib import Path

    output = Path("busco_sequences.tar.gz")

    archives = [
        Path("${multi_copy}"),
        Path("${single_copy}"),
        Path("${fragmented}"),
    ]

    # Write to a gz file as destination
    with tarfile.open(output, "w:gz") as dst:
        # For archive INSIDE the larger archive
        for archive in archives:
            # Read the gzip file
            with tarfile.open(archive, "r:gz") as src:
                for member in src:
                    # Rename the member to a child of the busco_sequences dir
                    # This is to satisfy the blobtoolkit tool
                    member.name = f"busco_sequences/{member.name}"
                    dst.addfile(member, src.extractfile(member))
    PY """
}
