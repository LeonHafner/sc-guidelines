
def addRuntimeToMeta(meta, sc_anndata, pb_anndata, time_file, method) {
    elapsed_time = time_file.text.readLines().head().trim()
    meta = meta + [(method): elapsed_time]
    return [meta, sc_anndata, pb_anndata]
}

def addRuntimeToMeta(meta, anndata, time_file, method) {
        elapsed_time = time_file.text.readLines().head().trim()
        meta = meta + [(method): elapsed_time]
        return [meta, anndata]
}