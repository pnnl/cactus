# How to contribute

Hello there. Looks like you are interested in contributing to [cactus](https://github.com/pnnl/cactus), which is just great, so welcome!
  
We are open to any kind of contribution: code, documentation, bug-fixing (even just reports), feature-requests and anything else you may come up with.

## Getting started

### Discussion

Discussions should be the first line of contact for new users to provide direct feedback on the library.
Post your Ideas for new features or examples, showcase your work using cactus in Show and tell, or get support for usage in Q&As, please post them in one of our categories.

## Contributing

### Code

To contribute code to the project:

1. Fork this repo
2. Set up the coding environment
3. Commit your code
4. Submit a pull request. It will be reviewed by maintainers and they'll give you proper feedback so you can iterate over it.

#### Considerations
- Make sure existing tests pass
- Make sure your new code is properly tested and fully-covered
- Following [The seven rules of a great Git commit message](https://chris.beams.io/posts/git-commit/#seven-rules) is highly encouraged
- When adding a new feature, branch from [main-branch](<project-main-branch>)


### Testing

As mentioned above, existing tests must pass and new features are required to be tested and fully-covered.

### Documenting

Code should be self-documented. But, in case there is any code that may be hard to understand, it must include some comments to make it easier to review and maintain later on.

## Development Environment

To make a contribution, you will need to set up the development environment. This repo uses [Rye](https://github.com/astral-sh/rye), but any `virtualenv` python tool will work (Pyenv, virtualenv, venv, etc.). 

### Using Rye ✅

> Make sure to have `rye` installed before performing these steps.

If using rye, you can quickly spin up the development environment using the following:

```shell
git clone git+https://github.com/pnnl/cactus.git
cd cactus
rye sync
```
> Note that you must first uninstall any current version you may already have.

This will will install the code in an editable configuration and install all the required dev dependencies.

### Not using Rye ❌

First you must install the code via:
```shell
git clone git+https://github.com/pnnl/cactus.git
cd cactus
pip install -e .
```
> Note that you must first uninstall any current version you may already have.

Then you can make sure you have all the dev dependencies installed by running
```shell
pip install -r requirements-dev.lock
```

## Linting, Formatting, and Testing 🧹🎨🧪

When making commits to your branch to prepare for your PR, make sure to follow proper code quality assessments by frequently running the linter, formatter, and test suite on your code. The PR will not be able to be merged if the linter and tests do not pass. 

> Be in the root directory of the repo before running these instructions.

### Using Rye ✅

To lint your code, use the command:
```shell
rye lint
rye run mypy --ignore-missing-imports .
```
This will return a list of current issues with your code.

---
To format your code, use the command:
```shell
rye fmt
```

---
To run the test suite on your code, use the command:
```bash
rye test
```

### Not using Rye ❌

To lint your code, use the command:
```shell
ruff check .
mypy --ignore-missing-imports .
```
This will return a list of current issues with your code.

---
To format your code, use the command:
```shell
ruff format .
```

---
To test your code, use the command:
```shell
pytest
```


