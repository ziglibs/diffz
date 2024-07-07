import os
import git

# Variables
REPO_PATH = '/Users/atman/code/opp/ziglibs/diffz'
FILE_NAME = 'DiffMatchPatch.zig'
OUTPUT_DIR = 'file-versions'

# Initialize the repository
repo = git.Repo(REPO_PATH)

# Create the output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Get a list of all commits that modified the file
commits = list(repo.iter_commits(paths=FILE_NAME))
commits.reverse()

# Loop through each commit
for i, commit in enumerate(commits):
    # Checkout the file from the specific commit
    file_content = (repo.git.show(f'{commit.hexsha}:{FILE_NAME}'))
    # Write the file content to the output directory with a suffix
    with open(os.path.join(OUTPUT_DIR, f'file-{i+1:02d}.zig'), 'w') as f:
        f.write(file_content)