Basic documentation to setup and use Meta's LLAMA via wexac resources;


The below is more or less copy-paste from https://aituts.com/llama/
## Step 1
```
conda create -n textgen python=3.10.9
conda activate textgen
pip install torch==1.13.1+cu116 torchvision==0.14.1+cu116 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu116
```

## Step 2
```
git clone https://github.com/oobabooga/text-generation-webui.git
cd text-generation-webui
pip install -r requirements.txt
pip install torch==1.12+cu113 -f https://download.pytorch.org/whl/torch_stable.html
mkdir repositories
cd repositories
git clone https://github.com/qwopqwop200/GPTQ-for-LLaMa.git
cd GPTQ-for-LLaMa
pip install ninja
conda install -c conda-forge cudatoolkit-dev
python setup_cuda.py install
```

## Step 3
```
pip install quant_cuda-0.0.0-cp310-cp310-win_amd64.whl
```

## Run as webapp
```
python server.py --cai-chat --model llama-7b --no-stream
```


