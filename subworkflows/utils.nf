def addRuntimeToMeta(meta, time_file, method) {
        elapsed_time = time_file.text.readLines().head().trim()
        meta = meta + [(method): elapsed_time]
        return [meta]
}

def addRuntimeToMeta(meta, anndata, time_file, method) {
        elapsed_time = time_file.text.readLines().head().trim()
        meta = meta + [(method): elapsed_time]
        return [meta, anndata]
}