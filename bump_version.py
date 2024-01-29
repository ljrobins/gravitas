import toml

def bump_version_pyproject_toml(file_path):
    # Read the TOML file
    with open(file_path, 'r') as file:
        data = toml.load(file)

    # Extract the current version
    current_version = data['project']['version']
    major, minor, patch = map(int, current_version.split('.'))

    # Increment the patch version
    new_version = f'{major}.{minor}.{patch+1}'

    # Update the version in the dictionary
    data['project']['version'] = new_version

    # Write the updated data back to the TOML file
    with open(file_path, 'w') as file:
        toml.dump(data, file)

    print(f"Version updated to {new_version}")

# Usage
bump_version_pyproject_toml('pyproject.toml')
