process PREPARE_FILES {

    container 'leonhafner/python'

    input:
    path lung_cancer_atlas

    output:
    path "splitted/*"

    script:
    """
    unzip ${lung_cancer_atlas} -d splitted
    """

    stub:
    """
    mkdir splitted
    touch splitted/1_file.h5ad
    touch splitted/2_file.h5ad
    touch splitted/3_file.h5ad
    """
}