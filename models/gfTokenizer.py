from geneformer import EmbExtractor

from geneformer.in_silico_perturber import load_model

if __name__ == "__main__":
    embex = EmbExtractor(model_type="Pretrained",
                         max_ncells=2000,
                         emb_layer=0,
                         emb_label=["region","cell_type"],
                         labels_to_plot=["cell_type"],
                         forward_batch_size=200,
                         nproc=4)

    model_dir = "/home/comp/cszrwang/project/Geneformer_test/Geneformer_model/geneformer-12L-30M"

    model = load_model(model_type="Pretrained", num_classes=0, model_directory=model_dir)

    print(model)

    # embs = embex.extract_embs("/home/comp/cszrwang/project/Geneformer_test/Geneformer_model/geneformer-12L-30M",
    #                           "/home/comp/cszrwang/project/Geneformer_test/STARmap/20180417_BZ5_control.dataset",
    #                           "/home/comp/cszrwang/project/Geneformer_test/STARmap/20180417_BZ5_control_output",
    #                           "output_prefix")



