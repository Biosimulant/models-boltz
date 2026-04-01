# biosim-boltz-pack-template

A template for creating additional Boltz-family module packs in this
repository.

## Installation

```bash
pip install -e .
```

## Usage

### In Python

```python
from src.my_pack.modules import ExampleBoltzModule

module = ExampleBoltzModule()
module.advance_to(0.01)
print(module.get_outputs())
```

### In YAML Space References

```yaml
models:
  - repo: Biosimulant/models-boltz
    alias: example_boltz
    manifest_path: models/example-boltz-module/model.yaml
    parameters: {}
```

## Development

```bash
pip install -e ".[dev]"
pytest
```

## Creating Your Own Pack

1. Copy this template.
2. Rename `my_pack` and the example module.
3. Update the entrypoint and metadata.
4. Keep the public module interface aligned with repo examples.
5. Add subprocess-boundary tests if the module shells out.
