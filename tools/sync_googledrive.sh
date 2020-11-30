#!/usr/bin/env bash

check_update_py=/Volumes/Huitian/Projects/CD8_DEV_SC/tools/check_up_to_date.py

git_repo_dir=/Volumes/Huitian/Projects/CD8_DEV_SC
google_dir=/Volumes/Huitian/GSuite_Scripps/x_Active_projects/CD8_DEV_SC

git_repo_figures=${git_repo_dir}/9_Figures
git_repo_updateFigures=${git_repo_dir}/z_update_figures
git_repo_updateManuscript=${git_repo_dir}/z_update_manuscript
git_repo_updateSlides=${git_repo_dir}/z_update_slides

google_figures=${google_dir}/9_Figures
google_updateFigures=${google_dir}/z_update_figures
google_updateManuscript=${google_dir}/z_update_manuscript
google_updateSlides=${google_dir}/z_update_slides

cd $git_repo_figures
cp -r -n * $google_figures

cd $git_repo_updateFigures
cp -r -n * $google_updateFigures

cd $git_repo_updateManuscript
cp -r -n * $google_updateManuscript

cd $git_repo_updateSlides
cp -r * $google_updateSlides

python $check_update_py $git_repo_updateFigures
python $check_update_py $google_updateFigures
