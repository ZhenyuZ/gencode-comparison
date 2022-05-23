VERSION = 1.0.0
REPO = bio-template
BRANCH_NAME?=unknown

init: init-hooks init-secrets

init-hooks:
	@echo
	@echo -- Installing Precommit Hooks --
	pre-commit install

init-secrets:
	@echo
	detect-secrets scan --update .secrets.baseline
	detect-secrets audit .secrets.baseline
