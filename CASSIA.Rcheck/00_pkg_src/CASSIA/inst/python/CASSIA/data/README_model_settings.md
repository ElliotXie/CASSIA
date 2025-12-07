# Model Settings Configuration

## Location

The `model_settings.json` file is located in the CASSIA Python package data directory:
```
CASSIA_python/CASSIA/data/model_settings.json
```

## Why This Location?

1. **Package Distribution**: The file is included with the package installation
2. **Version Control**: Updates to model settings are versioned with the package
3. **Accessibility**: Both Python and R implementations can access it
4. **Consistency**: Ensures all users have the same model configurations

## How It Works

The model settings system looks for the JSON file in this order:

1. **Package data directory** (preferred): `CASSIA/data/model_settings.json`
2. Same directory as `model_settings.py`
3. Parent directories (for development)
4. Project root
5. **Fallback**: Uses built-in configuration if no file found

## Customization

If you want to customize model settings:

1. **For development**: Place your custom `model_settings.json` in the project root
2. **For production**: Modify the file in the package data directory
3. **For temporary use**: Pass a custom path to the `ModelSettings` constructor

## Benefits

- ✅ Always available with package installation
- ✅ Consistent across all users
- ✅ Version controlled with the package
- ✅ Easy to update in future releases
- ✅ Accessible from both Python and R implementations