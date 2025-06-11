#!/bin/bash

# Makefile for concord project
# $0 compile  / c      Compile the project
# $0 build    / b      Build the project
# $0 run      / r      Run the project
# $0 test     / t      Run tests
# $0 docs     / d      Compile mdbook documentation
# $0 release [type]    Mark as releaser and create a new release

project_name=$(grep -Po 'set\s*\(\s*project_name\s+\K[^)]+' CMakeLists.txt)
project_cap=$(echo "$project_name" | tr '[:lower:]' '[:upper:]')

if [[ -z "$project_name" ]]; then
    echo "Error: project_name not found in CMakeLists.txt"
    exit 1
fi

echo -e "------------------------------------------"
echo -e "Project: $project_name"
echo -e "------------------------------------------\n"


# @cmd compile project
# @alias c
compile() {
    CURR_DIR=$(pwd)
    if [[ ! -d "$TOP_HEAD/build" ]] then
        mkdir "$TOP_HEAD/build";
    else
        rm -rf "$TOP_HEAD/build"
        mkdir "$TOP_HEAD/build"
    fi
    cd "$TOP_HEAD/build"
    echo "cmake -Wno-dev -D${project_cap}_BUILD_EXAMPLES=ON -D${project_cap}_ENABLE_TESTS=ON .."
    cmake -Wno-dev -D${project_cap}_BUILD_EXAMPLES=ON -D${project_cap}_ENABLE_TESTS=ON ..
    cd "$CURR_DIR"
}


# @cmd build project
# @alias b
build() {
    CURR_DIR=$(pwd)
    cd "$TOP_HEAD/build"
    make -j$(nproc) || true
    cd "$CURR_DIR"
}


# @cmd run project
# @alias r
run() {
    CURR_DIR=$(pwd)
    $TOP_HEAD/build/./main
    cd "$CURR_DIR"
}

# @cmd run tests
# @alias t
test() {
    CURR_DIR=$(pwd)
    cd "$TOP_HEAD/build"
    ctest --verbose --output-on-failure || true
    cd "$CURR_DIR"
}

# @cmd compile mdbook
# @arg type![mdbook|doxygen] Documentation type
docs() {
    case $argc_type in
        mdbook)
            if ! command -v mdbook &> /dev/null; then
                echo "mdbook is not installed. Please install it first."
                exit 1
            else
                mdbook build $TOP_HEAD/book --dest-dir $TOP_HEAD/docs
                git add --all && git commit -m "docs: building website/mdbook"
            fi
            ;;
        doxygen)
            if ! command -v doxygen &> /dev/null; then
                echo "doxygen is not installed. Please install it first."
                exit 1
            fi
            ;;
        *)
            echo "Invalid documentation type. Use 'mdbook' or 'doxygen'."
            exit 1
            ;;
    esac
}

# @cmd release project
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



eval "$(argc --argc-eval "$0" "$@")"
