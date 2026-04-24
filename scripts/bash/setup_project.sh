# Ask for project name
read -p "Enter project name: " project_name

# Create project structure
mkdir -p \
"$project_name"/data/raw \
"$project_name"/data/processed \
"$project_name"/output \
"$project_name"/scripts/bash \
"$project_name"/scripts/qmd

touch "$project_name/$project_name.Rproj"
# Confirmation message
echo "Project '$project_name' created successfully."
