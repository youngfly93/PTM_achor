# Makefile for PTM-HLA Anchor Coupling Analysis Pipeline

# 默认目标
.PHONY: all
all: meta

# 环境设置
.PHONY: setup
setup:
	bash setup_env.sh

# 构建元数据
.PHONY: meta
meta:
	python3 build_meta.py
	python3 update_hla_config.py

# 测试运行
.PHONY: test
test: meta
	python3 run_pipeline.py --meta all_meta_test.tsv --test

# 批次分析测试
.PHONY: test_batch
test_batch: meta
	python3 run_pipeline_batch.py --meta all_meta_test.tsv --test --mode batch

# 传统分析（所有数据合并）
.PHONY: traditional
traditional: meta
	python3 run_pipeline.py --meta all_meta_updated.tsv

# 批次分析（推荐）
.PHONY: batch
batch: meta
	python3 run_pipeline_batch.py --meta all_meta_updated.tsv --mode batch

# 批次分析 + MHCflurry
.PHONY: batch_mhc
batch_mhc: meta
	python3 run_pipeline_batch.py --meta all_meta_updated.tsv --mode batch --mhcflurry

# 清理结果
.PHONY: clean
clean:
	rm -rf results/ results_batch/ results_hla_test/ results_full/
	rm -f nohup.out

# 清理所有（包括元数据）
.PHONY: clean_all
clean_all: clean
	rm -f all_meta.tsv all_meta_updated.tsv all_meta_test.tsv
	rm -f hla_reference_stats.txt

# 生成报告
.PHONY: report
report:
	@echo "=== Pipeline Results Summary ==="
	@if [ -f results_batch/pipeline_summary.txt ]; then \
		cat results_batch/pipeline_summary.txt; \
	elif [ -f results/pipeline_summary.txt ]; then \
		cat results/pipeline_summary.txt; \
	else \
		echo "No results found. Run 'make batch' or 'make traditional' first."; \
	fi

# 检查异质性
.PHONY: heterogeneity
heterogeneity:
	@if [ -f results_batch/heterogeneity_report.txt ]; then \
		cat results_batch/heterogeneity_report.txt; \
	else \
		echo "No heterogeneity report found. Run 'make batch' first."; \
	fi

# 帮助信息
.PHONY: help
help:
	@echo "PTM-HLA Anchor Coupling Analysis Pipeline"
	@echo "========================================="
	@echo ""
	@echo "Setup commands:"
	@echo "  make setup         - Set up conda environment"
	@echo "  make meta          - Build metadata files"
	@echo ""
	@echo "Analysis commands:"
	@echo "  make test          - Run test analysis (5 samples)"
	@echo "  make test_batch    - Run batch test analysis"
	@echo "  make traditional   - Run traditional analysis (all data merged)"
	@echo "  make batch         - Run batch analysis with meta-analysis (recommended)"
	@echo "  make batch_mhc     - Run batch analysis with MHCflurry"
	@echo ""
	@echo "Utility commands:"
	@echo "  make report        - Show analysis summary"
	@echo "  make heterogeneity - Show heterogeneity assessment"
	@echo "  make clean         - Remove result files"
	@echo "  make clean_all     - Remove all generated files"
	@echo "  make help          - Show this help message"
	@echo ""
	@echo "Quick start:"
	@echo "  make setup && make batch"