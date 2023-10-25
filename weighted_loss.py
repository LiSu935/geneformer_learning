import numpy as np
import torch
from torch import nn

adata = None
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

celltype_id_labels = adata.obs["celltype"].astype("category").cat.codes.values
class_num = np.unique(celltype_id_labels, return_counts=True)[1].tolist()
class_weight = torch.tensor([(1 - (x / sum(class_num))) ** 2 for x in class_num])

criterion_cls = nn.CrossEntropyLoss(weight=class_weight.to(device))
