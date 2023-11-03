import numpy as np
import torch
import torch.nn as nn
from transformers import Trainer

# load the ms adata
adata = None
celltype_key = 'celltype'

# calculate weights
celltype_id_labels = adata.obs[celltype_key].astype("category").cat.codes.values
class_num = np.unique(celltype_id_labels, return_counts=True)[1].tolist()
class_weight = torch.tensor([(1 - (x / sum(class_num))) ** 2 for x in class_num])

# Add weight to the cross entropy loss function
weighted_loss = nn.CrossEntropyLoss(weight=class_weight)


class MyTrainer(Trainer):
    # overwrite the function of compute_loss
    def compute_loss(self, model, inputs, return_outputs=False):
        outputs = model(**inputs)

        # put loss to device
        loss_fct = weighted_loss.to(outputs['logits'].device)

        # calculate weighted loss
        loss = loss_fct(outputs['logits'], inputs['labels'])

        return (loss, outputs) if return_outputs else loss
