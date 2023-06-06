import os
import pandas as pd
import dendropy
from tqdm import tqdm
from dendropy.calculate import treecompare

if __name__ == "__main__":
    d1 = "./test/trees1/"
    d2 = "./test/trees2/"

    res = []
    for file in tqdm(os.listdir(d1)):
        tns = dendropy.TaxonNamespace()
        tree1 = dendropy.Tree.get_from_path(
            os.path.join(d1, file),
            "newick",
            taxon_namespace=tns,
        )
        tree2 = dendropy.Tree.get_from_path(
            os.path.join(d2, file),
            "newick",
            taxon_namespace=tns,
        )

        tree1.encode_bipartitions()
        tree2.encode_bipartitions()

        res.append(
            {
                "id": file.split(".")[0],
                "RF": treecompare.symmetric_difference(tree1, tree2),
                "wRF": treecompare.weighted_robinson_foulds_distance(tree1, tree2),
            }
        )

    df = pd.DataFrame(res).set_index("id")
    df.to_csv("python.dp.csv")
