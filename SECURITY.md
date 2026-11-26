# Security policy

## Scope

PyYAeHMOP is research software that processes user-provided geometry and parameter files, runs numerical Python computations, and optionally launches a desktop GUI.

## Supported version

The currently documented source version is **0.4.0**. Security fixes, if required, will be handled in subsequent source updates. There is no guaranteed commercial support or response timeline.

## Reporting a vulnerability

**Do not publish exploitable security details or private data in a public issue.** If GitHub Private Vulnerability Reporting is available for this repository, use the **Security → Report a vulnerability** interface. Otherwise, contact the repository maintainer privately through the GitHub account profile to coordinate disclosure.

Please include affected version/commit, platform, minimal reproduction steps and impact. Do not send credentials or sensitive structures.

## Operating safely

- Treat untrusted structure files and parameter paths as inputs that require review.
- The GUI starts a child Python process to perform calculations; use it with input sources you trust.
- Use a dedicated Python environment for scientific dependencies.
- Do not run contributed shell scripts with root permissions without review.
- Avoid processing confidential structures in shared result folders.
- GitHub workflow files in this repository use manual dispatch only; no computation is automatically triggered on push.
