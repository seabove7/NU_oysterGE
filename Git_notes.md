## Working on GitHub Guide

**Suggestions compiled by [Colleen Bove](http://colleenbove.science) for summer 2022 BRITE intern collaborations**

<br/>


### Getting Started

1) Open up your project in RStudio and **ALWAYS** start by checking the status of the project `git status`.

2) If your branch is up to date with origin/main, then you can go ahead and pull (`git pull`) from main and start working! If your branch is behind the origin/main, you will need to remedy these changes before pulling.
  
   - You hay have a `.DS_Store` file with changes that have not been committed to main. I generally just restore these from the origin/main branch like so: `git restore <filename>`. This will overwrite what you have locally with the current origin/main file so be careful with this, but it generally safe with `.DS_Store` files. 
   
   - Another issue is if you have been editing any repository files on your local computer before pulling. If this is the case, you will need to work to merge these so no edits are lost from either end. This is complicated and will take me than this tutorial will cover right now.   

3) Work on editing your code and any other files in your repo! If you are doing a lot of major edits, I **highly** recommend making several commits (`git commit`) throughout the process to track these major changes. For example, I will generally commit changes after added a new figure or new analysis.    

4) When you are finished working on that project for a while or especially for the day, you should commit and push **ALL** changes to git (see commands below, but you will follow this generally: `git add *`, `git commit -m "COMMENT"`, `git push`).

5) After you push all changes, I recommend you check your [GitHub](https://github.com) repo online to make sure the push went through successfully. If it did, close all files in RStudio and close the project to make sure you do not accidentally make changes to the files when you do not want to. Then you are all set!

<br/>


### Some helpful commands:  

1. Check the status of your git repo with main  
`git status`

2. Pull the most up-to-date edits from the main branch (and ensure you are fully up to date)  
`git pull`

3. Add all updates to the code (you can either specify the file or use * to add all changed files)  
`git add <filename>`
`git add *`

4. Commit all staged edits with a short comment on what was done (e.g., "added figure for PCAs and PERMANOVA")  
`git commit -m "COMMENT ON UPDATE"`

5. Push the committed changes to the main branch  
`git push`

6. If you have files that have been modified that are preventing a pull request (and that are not as important i.e., the .DS_Store files), you can use the following:  
`git restore <filename>`

<br/>


**This is just a basic intro to integrating git with Rmarkdown for version control and better collaboration on coding. There are many 'Best Practices' guides available across the internet that are also helpful with this! [Here is one from GitHub with some basic information](https://docs.github.com/en/get-started/using-git/about-git)**

<br/>
<br/>
