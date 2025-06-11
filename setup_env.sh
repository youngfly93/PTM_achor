#!/bin/bash
# 修饰位点 × HLA 锚定位点耦合分析环境设置

echo "Creating conda environment for anchor-modification analysis..."

# 创建conda环境
conda create -n anchor_mod python=3.10 -y

# 激活环境
source $(conda info --base)/etc/profile.d/conda.sh
conda activate anchor_mod

# 安装基础包
conda install -y pandas pyteomics tqdm seaborn matplotlib scipy

# 安装mhcflurry (HLA结合预测)
pip install mhcflurry logomaker

# 下载mhcflurry模型 (首次使用需要)
python -c "
from mhcflurry import Class1PresentationPredictor
print('Downloading MHCflurry models...')
predictor = Class1PresentationPredictor.load()
print('MHCflurry models downloaded successfully!')
"

echo "Environment setup completed!"
echo "To activate: conda activate anchor_mod"

# 可选：如果有netMHCpan，检查PATH
if command -v netMHCpan &> /dev/null; then
    echo "netMHCpan found in PATH: $(which netMHCpan)"
else
    echo "netMHCpan not found in PATH - using mhcflurry only"
fi