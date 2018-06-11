
cp .gitmodules .gitmodules.old

sed -e 's/git@github\.com\:sim82/git\:\/\/github\.com\/sim82/g' .gitmodules.old > .gitmodules
