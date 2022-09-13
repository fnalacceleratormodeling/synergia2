# Useful Github Stuff

## Debugging CI workflows

One can access the CI runners with [tmate](https://github.com/marketplace/actions/debugging-with-tmate). These hooks can be enabled when manually launching a workflow that allow one to connect to the GithubActions runner at predefined checkpoints. At the moment, the following workflows have this feature enabled:
- ubuntu-clang
- macos-clang

## Checking workflow action `yml` files

Static [github action syntax checker](https://rhysd.github.io/actionlint)
