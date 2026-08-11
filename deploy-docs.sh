#!/bin/bash
# ====================================================================
# Script to build IDLdoc locally and push directly from out-docs/
# ====================================================================

# 1. Clean out any old build files from out-docs to prevent duplication
rm -rf out-docs
mkdir -p out-docs

# 2. Run local IDLdoc by piping commands into IDL.
# (Now out-docs/ is clean, and since it is empty, IDLdoc won't find old files)
idl << 'EOF'
idldoc, ROOT='.', OUTPUT='out-docs/', TITLE='SSTRED IDL Code Documentation'
EXIT
EOF

# 3. Safety check: Ensure the output folder actually contains data
if [ ! -d "out-docs" ] || [ -z "$(ls -A out-docs)" ]; then
    echo "Error: Folder out-docs/ was not created or is empty. Aborting."
    exit 1
fi

# 4. Jump directly into out-docs/ using pushd
pushd out-docs > /dev/null

# 5. Initialize a clean, isolated Git space inside the output folder
git init
git checkout -b gh-pages

# 6. Commit the fresh HTML files locally
git add .
git commit -m "Automated local IDLdoc update"

# 7. Force-push the clean HTML content to the hidden branch on GitHub
git push -f git@github.com:ISP-SST/sstred.git gh-pages

# 8. Jump back to your project root folder
popd > /dev/null

# 9. Clean up out-docs/ completely so it doesn't take up space or interfere locally
rm -rf out-docs

echo ""
echo "Success! Your SSTRED documentation has been thoroughly cleaned, rebuilt,"
echo "and pushed to GitHub Pages without any duplicate files."
echo "You are safely back in your project directory: $(pwd)"
