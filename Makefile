# General
PACKAGE_NAME=ACTIONet
CURR_VERSION=$(shell head -n1 ACTIONet/__init__.py | cut -f3 -d' '| sed 's/"//g')
SOURCE_CONDA=source $$(conda info --base)/etc/profile.d/conda.sh

# Testing
COVERAGE_REPORT_FILE=coverage.xml
PYTEST_REPORT_FILE=report.xml

## Linting
.PHONY: install_pre_commit
install_pre_commit:
	pre-commit install

.PHONY: format
format: install_pre_commit
	pre-commit run black --all-files
	pre-commit run isort --all-files

.PHONY: type_check
type_check:	install_pre_commit
	pre-commit run mypy --all-files

.PHONY: lint
lint:	install_pre_commit
	pre-commit run --all-files

## Installation
.PHONY: update_env install_env install install_dev reinstall resintall_dev
update_env:
	$(SOURCE_CONDA) && conda env update -f environment.yaml || conda env create -f environment.yaml

install_env:
	$(SOURCE_CONDA) && conda deactivate && conda env remove -n actionet
	$(SOURCE_CONDA) && conda deactivate && conda env create -f environment.yaml

install: install_env
	$(SOURCE_CONDA) && conda activate actionet && pip install --find-links https://conda.aws.insitro.com/pypi/ actionet==$(CURR_VERSION)

install_dev:	install_env
	$(SOURCE_CONDA) && conda activate actionet && pip install --find-links https://conda.aws.insitro.com/pypi/ -e .

reinstall:	update_env
	$(SOURCE_CONDA) && conda activate actionet && pip install --find-links https://conda.aws.insitro.com/pypi/ actionet==$(CURR_VERSION)

reinstall_dev:	update_env
	$(SOURCE_CONDA) && conda activate actionet && pip install --find-links https://conda.aws.insitro.com/pypi/ -e .

## Testing
.PHONY: pytest
pytest:
	pytest -v --junitxml $(PYTEST_REPORT_FILE) $(PACKAGE_NAME)

.PHONY: pytest-cov
pytest-cov:
	pytest --cov $(PACKAGE_NAME) $(PACKAGE_NAME)

.PHONY: pytest-integration
pytest-integration:
	pytest $PACKAGE_NAME -m integration_test

## Publish docs
.PHONY: build_docs
build_docs:
	cd docs && rm -r build
	cd changelog && chmod +x *sh && bash ./generate_change_log.sh
	cd docs && make html

