#!/bin/bash


# @cmd build cmake
# @alias b
# # @flag      --arrow        Enable gui support
buildit() {
    CURR_DIR=$(pwd)
    if [[ ! -d "$TOP_HEAD/build" ]] then
        mkdir "$TOP_HEAD/build";
    else
        rm -rf "$TOP_HEAD/build"
        mkdir "$TOP_HEAD/build"
    fi
    cd "$TOP_HEAD/build"
    if ! [ -z "$argc_arrow" ]; then
        cmake -Wno-dev ..
    else
        cmake -Wno-dev -DRERUN_DOWNLOAD_AND_BUILD_ARROW=OFF ..
    fi
    cd "$CURR_DIR"
}


# @cmd make project
# @alias m
makeit() {
    CURR_DIR=$(pwd)
    cd "$TOP_HEAD/build"
    make
    cd "$CURR_DIR"
}


# @cmd mark as releaser
# @alias r
# @arg type![patch|minor|major] Release type
release() {
    CURRENT_VERSION=$(grep -E '^project\(.*VERSION [0-9]+\.[0-9]+\.[0-9]+' CMakeLists.txt \
        | sed -E 's/.*VERSION ([0-9]+\.[0-9]+\.[0-9]+).*/\1/')
    IFS='.' read -r MAJOR MINOR PATCH <<< "$CURRENT_VERSION"
    case $argc_type in
        major)
            MAJOR=$((MAJOR + 1))
            MINOR=0
            PATCH=0
            ;;
        minor)
            MINOR=$((MINOR + 1))
            PATCH=0
            ;;
        patch)
            PATCH=$((PATCH + 1))
            ;;
    esac
    version="$MAJOR.$MINOR.$PATCH"
    sed -i -E "s/(project\(.*VERSION )[0-9]+\.[0-9]+\.[0-9]+/\1$version/" CMakeLists.txt
    git cliff --tag $version > CHANGELOG.md
    changelog=$(git cliff --unreleased --strip all)
    git add -A && git commit -m "chore(release): prepare for $version"
    echo "$changelog"
    git tag -a $version -m "$version" -m "$changelog"
    git push --follow-tags --force --set-upstream origin main
    gh release create $version --notes "$changelog"
}


# @cmd compile mdbook
# @option    --dest_dir <dir>    Destination directory
mdbook() {
    mdbook build $TOP_HEAD/book --dest-dir $TOP_HEAD/docs
    git add --all && git commit -m "docs: building website/mdbook"
}


eval "$(argc --argc-eval "$0" "$@")"
