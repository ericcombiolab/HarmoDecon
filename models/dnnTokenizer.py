import torch
import torch.nn as nn

class DnnTokenizer(nn.Module):
    def __init__(self, input_dimension):
        super(DnnTokenizer, self).__init__()

        # Define the layers of your DNN
        self.fc1 = nn.Linear(input_dimension, 1024)
        self.relu1 = nn.ReLU()
        self.fc2 = nn.Linear(1024, 512)
        self.relu2 = nn.ReLU()
        self.LayerNorm = nn.LayerNorm(512, eps=1e-12)
        self.fc3 = nn.Linear(512, 512)
        self.relu3 = nn.ReLU()
        self.LayerNorm = nn.LayerNorm(512, eps=1e-12)

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu1(x)
        x = self.fc2(x)
        x = self.relu2(x)
        x = self.fc3(x)
        x = self.relu3(x)
        x = self.LayerNorm(x)
        return x

if __name__ =="__main__":
    x_genes = torch.rand((1, 200,20000))
    tokenizer = DnnTokenizer(input_dimension=x_genes.shape[-1])
    x_token = tokenizer(x_genes)
    print(x_token.shape)