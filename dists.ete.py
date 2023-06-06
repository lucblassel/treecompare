import os
import pandas as pd
from tqdm import tqdm
from ete3 import Tree

if __name__ == "__main__":
    d1 = "./test/trees1/"
    d2 = "./test/trees2/"

    res = []
    for file in tqdm(os.listdir(d1)):
        tree1 = Tree(
            os.path.join(d1, file),
        )
        tree2 = Tree(
            os.path.join(d2, file),
        )

        compare = tree1.compare(tree2)

        res.append(
            {
                "id": file.split(".")[0],
                "RF": compare["rf"],
                "nRF": compare["norm_rf"],
            }
        )

    df = pd.DataFrame(res).set_index("id")
    df.to_csv("python.dp.csv")
