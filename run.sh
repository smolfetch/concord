#!/bin/bash

# @cmd compile project
# @alias c
compileit() {
    CURR_DIR=$(pwd)
    if [[ ! -d "$TOP_HEAD/build" ]] then
        mkdir "$TOP_HEAD/build";
    else
        rm -rf "$TOP_HEAD/build"
        mkdir "$TOP_HEAD/build"
    fi
    cd "$TOP_HEAD/build"
    cmake -Wno-dev -DBUILD_EXAMPLES=ON ..
    cd "$CURR_DIR"
}


# @cmd build project
# @alias b
buildit() {
    CURR_DIR=$(pwd)
    cd "$TOP_HEAD/build"
    make
    cd "$CURR_DIR"
}


# @cmd run project
# @alias r
runit() {
    CURR_DIR=$(pwd)
    $TOP_HEAD/build/./main
    cd "$CURR_DIR"
}

# @cmd mark as releaser
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
    if [ -n "$LATEST_TAG" ]; then
        # Get changelog content for release notes (changes since last tag)
        changelog=$(git cliff $LATEST_TAG..HEAD --strip all)
        # Generate changelog and prepend to existing file (changes since last tag)
        git cliff --tag $version $LATEST_TAG..HEAD --prepend CHANGELOG.md
    else
        # First release - get all changes
        changelog=$(git cliff --unreleased --strip all)
        git cliff --tag $version --unreleased --prepend CHANGELOG.md
    fi
    sed -i -E "s/(project\(.*VERSION )[0-9]+\.[0-9]+\.[0-9]+/\1$version/" CMakeLists.txt
    git add -A && git commit -m "chore(release): prepare for $version"
    echo "$changelog"
    git tag -a $version -m "$version" -m "$changelog"
    git push --follow-tags --force --set-upstream origin develop
    gh release create $version --notes "$changelog"
}


# @cmd compile mdbook
# @option    --dest_dir <dir>    Destination directory
mdbook() {
    mdbook build $TOP_HEAD/book --dest-dir $TOP_HEAD/docs
    git add --all && git commit -m "docs: building website/mdbook"
}


eval "$(argc --argc-eval "$0" "$@")"
