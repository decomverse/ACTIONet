git submodule add --branch core --depth 1 -- https://github.com/shmohammadi86/ACTIONet.git src/ACTIONet
git config -f .gitmodules submodule.src/ACTIONet.shallow true
git submodule update --init --recursive --depth 1 --
git submodule update --remote --depth 1 --

