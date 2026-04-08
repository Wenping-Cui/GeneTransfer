# Repo Cleanup Prompts

A reusable set of prompts for structurising, testing, documenting, and maintaining a scientific code repository. Run them in order for a full cleanup, or pick individual steps as needed.

---

## 1. Structurize the codebase

```
Please structurize the codes as a professional coder. Keep the changes in a JSON file (config.json). 
Also provide test cases to run the code. Specifically:
- Separate configuration/default parameters into config.json
- Add clear docstrings and section headers to all functions
- Fix inconsistent indentation and formatting
- Fix any obvious bugs found during the review
- Create a test/ directory with a runtests file covering all exported functions
```

---

## 2. Detect and configure dev servers

```
Detect my project's dev servers and save all their configurations to .claude/launch.json, 
then ask which ones to start. Use this format:

{
  "version": "0.0.1",
  "configurations": [
    {
      "name": "<server-name>",
      "runtimeExecutable": "<command>",
      "runtimeArgs": ["<args>"],
      "port": <port>
    }
  ]
}

Use runtimeExecutable for the command (e.g. "yarn", "npm", "node", "python", "jupyter") 
and runtimeArgs for its arguments. Call preview_start for each server the user wants to run.
```

---

## 3. Check notebook consistency with source files

```
Go through all Jupyter notebooks in the repository and make sure they are consistent with 
the updated source files. Check for:
- Incorrect module names in import/using statements
- Outdated function signatures or renamed functions
- Functions exported from the module but called with the wrong name
- Any references to removed or renamed variables
Fix any inconsistencies found.
```

---

## 4. Verify refactored code produces identical output

```
Did you keep the previous version of the main script? 
If not, retrieve it from git history.
Write a comparison script that:
- Reproduces both the old and new versions of the core simulation/processing function
- Runs both with an identical random seed and realistic parameters matching the actual use case
- Compares outputs across all supported modes/configurations (e.g. constraint types)
- Reports whether outputs are bit-for-bit identical
Save the script to test/compare_versions.jl (or equivalent) so it can be re-run later.
```

---

## 5. Code review and cleanup (reuse, quality, efficiency)

```
Please review all changed files for reuse, quality, and efficiency, then fix any issues found.
Focus on:
- Duplicate functions or logic that can be unified
- Redundant allocations or computations (e.g. computing the same value twice)
- Missing bounds guards that could cause silent errors
- Dead variables or unused return values
- Functions that mix I/O with computation (separate them so the logic is testable)
- Near-identical code blocks that should share an abstraction
After reviewing, fix all issues and re-run the tests to confirm nothing broke.
```

---

## 6. Add a comprehensive README

```
Add a comprehensive README.md for this repository. The project aims for: <PAPER_URL>

Include:
- Full paper citation (title, authors, journal, DOI, year)
- Scientific overview: what problem does the paper solve, what is the key finding
- Repository structure: a file tree with one-line descriptions of each file/folder
- Model or algorithm description: key inputs, mechanisms, and outputs
- Requirements: language version and all dependencies
- Usage guide: a minimal working example with the most important CLI flags/parameters
- A parameter table for the most commonly adjusted settings
- Instructions for running tests
- Contact information
```

---

## 7. Commit changes

```
Commit all staged changes with a descriptive commit message that summarises:
- What was changed and why (not just what)
- List the key fixes, additions, and removals as bullet points
```

---

## 8. Review what changed since last commit

```
Can you take a look at what I did? 
Summarise all uncommitted changes across modified and untracked files, 
and tell me whether they look correct and consistent.
```

---

## 9. Full pipeline (run all steps at once)

Combine the above into a single prompt for a complete repo cleanup from scratch:

```
Please do a full repository cleanup:
1. Structurize all source files as a professional coder — consistent formatting, 
   docstrings, section headers, and bug fixes.
2. Extract all default parameters into config.json.
3. Write a comprehensive test suite in test/runtests.jl covering all exported functions.
4. Check all Jupyter notebooks for consistency with the updated source files and fix any issues.
5. Detect dev servers and save configurations to .claude/launch.json.
6. Review all changes for code reuse, quality, and efficiency — fix any issues found.
7. Verify the refactored code produces identical output to the original using a fixed random seed.
8. Write a comprehensive README.md. The paper this repo accompanies is: <PAPER_URL>
9. Commit everything with a descriptive commit message and push.
```

---

## Tips

- Always run the test suite after each step to catch regressions early.
- Pass a DOI or paper URL when requesting README generation — it lets Claude fetch the abstract automatically.
- For the version comparison script, use parameters that match the real use case (realistic N, T, L) rather than toy values — it makes the check more meaningful.
- The `config.json` approach works for any language: Julia (`JSON.jl`), Python (`json`), R (`jsonlite`), etc.
