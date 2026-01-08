# Regenerate all CYTools fixtures
# Requires: pip install cytools (or use docker)
# Run: make all

SCRIPTS_DIR := scripts
FIXTURES_DIR := fixtures

GENERATORS := \
	generate_glsm.py \
	generate_triangulation.py \
	generate_linear_relations.py \
	generate_divisor_basis.py \
	generate_intersection_full.py \
	generate_intersection_basis.py \
	generate_flat_direction.py \
	generate_filter_tensor.py

.PHONY: all clean

all: $(FIXTURES_DIR)
	@for script in $(GENERATORS); do \
		echo "Running $$script..."; \
		python3 $(SCRIPTS_DIR)/$$script || exit 1; \
	done
	@echo "All fixtures generated."

$(FIXTURES_DIR):
	mkdir -p $(FIXTURES_DIR)

clean:
	rm -f $(FIXTURES_DIR)/*.json
