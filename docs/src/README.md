# Getting started
Using Muscade, and developping new elements or materials for Muscade, requires writing and executing code in the Julia programming language. This can be done using high quality freeware. There are various options to do so, the recommended one is detailed below.  It is then necessary to install the Muscade package.

## Installing Julia

The following instructions have been tested on Windows:

First, from [julialang.org/downloads](https://julialang.org/downloads), download and install the "current stable release" of Julia. (no need to add Julia to path, but you can).

From [code.visualstudio.com](https://code.visualstudio.com), download and install the "stable build" of Visual Studio Code (VScode). Launch VCCode.  On the left bar , click the "extensions" icon (four squares, of which one is not touching the others).  On the top of the pane, in the search field, type `Julia`.  Install the "Julia (Julia language support)" extension (not the "Julia Insider"). The extension should find the Julia language (in Windows `C:\Users\philippem\AppData\Local\Programs\Julia-1.8.0\bin\julia.exe`).  If not, clicking on the installed Julia extension displays a "Details" test with instructions.

Create a new file, and save it with a title ending in `.jl`. If VSCode has detected the Julia language isntallation, a right pointing triangle allowing to execute the current file. In the file, type the canonical `print("Hello world!")`, and execute by clicking the triangle.  This should open a Julia REPL (a Julia terminal) pane at the bottom of the VSCode window.  After a little compiling, the greeting should appear (a second execution will be much faster, typical Julia).

Other extensions that are highly recommended are "Fast Unicode Math Character" (with it installed, typing `\alpha` and then pressing TAB will cause `α`"" to appear), `Julia Color Themes` for syntax highlighting, `vscode-js-profile-flame` (a part of producing graphs for profiling of code). GIT users should consider "GitLens".

## Installing Muscade

In the REPL, type `]` (to go into package management mode), followed by 

- `add Muscade`.
- Press the `backspace` key to leave the package manager.

