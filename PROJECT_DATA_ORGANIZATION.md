# Project Data Organization

This document describes the recommended folder structure and naming convention for organizing project data.

## Overview

Each project should have a clear, informative name that reflects the dataset or experiment. For example:

`hpc_ctx_project`

Within each project folder:

1. Create one subfolder per animal.
2. Within each animal folder, create one subfolder per recording session.
3. Within each session folder, store one subfolder per recording subepoch.
4. Within each subepoch folder, store the raw Intan data and associated video files.

## Folder Hierarchy

The expected hierarchy is:

```text
project_name/
  HP01/
    session_name/
      subepoch_name/
  HP02/
    session_name/
      subepoch_name/
```

Using the example project:

```text
hpc_ctx_project/
  HP18/
    hp18_day18_20250422/
      hp18_probe_250422_104827/
      hp18_presleep_250422_110045/
      hp18_cheeseboard1_250422_131312/
      hp18_postsleep1_250422_134641/
      hp18_cheeseboard2_250422_152803/
      hp18_postsleep2_250422_160304/
```

## Naming Rules

### 1. Project folder

Project folders should use informative names. Example:

`hpc_ctx_project`

This makes it easier to distinguish projects at a glance and reduces confusion when multiple datasets are stored on the same server.

### 2. Animal folder

Each animal gets its own folder, using a short ID such as:

`HP01`, `HP02`, `HP18`

Using one folder per animal keeps sessions grouped together and makes it easy to browse all data collected from that subject.

### 3. Session folder

Each recording session should have its own folder. Example:

`hp18_day18_20250422`

This session name contains:

- `hp18`: the animal ID
- `day18`: the recording day
- `20250422`: the calendar date in `YYYYMMDD` format

This naming format makes it easy to identify the animal, the experimental day, and the actual recording date directly from the folder name.

### 4. Subepoch folder

Within each session folder, create one subfolder for each recording epoch started and stopped by the experimenter. Examples:

- `hp18_probe_250422_104827`
- `hp18_presleep_250422_110045`
- `hp18_cheeseboard1_250422_131312`
- `hp18_postsleep1_250422_134641`
- `hp18_cheeseboard2_250422_152803`
- `hp18_postsleep2_250422_160304`

These names intentionally include:

- the animal ID
- the recording epoch label
- the date stamp
- the time stamp

The timestamp portion is generated automatically by the Intan software. Including the animal ID and day-specific context in the subepoch name helps prevent copy/paste mistakes when multiple animals are being recorded on the same day and data are transferred to the server later.

## What Goes Inside Each Subepoch Folder

Each subepoch folder should contain the raw acquisition outputs for that recording block, including:

- raw Intan data
- corresponding video data

The subepoch folder is the lowest-level organizational unit for raw data storage.

## Why This Structure Helps

This organization provides several benefits:

- it keeps data grouped first by project, then by animal, then by session, then by recording block
- it makes sessions easy to identify from folder names alone
- it reduces the chance of mixing data across animals or days
- it supports cleaner server transfers at the end of recording days
- it preserves a predictable structure for downstream analysis pipelines

## Example Full Path

An example session path is:

`U:\data\hpc_ctx_project\HP18\hp18_day18_20250422`

An example subepoch path inside that session is:

`U:\data\hpc_ctx_project\HP18\hp18_day18_20250422\hp18_probe_250422_104827`

