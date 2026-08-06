# runners

## Purpose

The `runners` package orchestrates Annotation QC workflows.

A runner connects the different components of the framework by executing them in the correct order.

Runners coordinate the workflow but should contain very little business logic.

---

# Responsibilities

A runner is responsible for

* selecting the appropriate parser
* loading annotation data
* calling one or more metric functions
* collecting metric outputs
* invoking the appropriate report generator
* handling workflow configuration

Think of a runner as the controller of the QC pipeline.

---

# What does NOT belong here?

Do not implement

* annotation parsing
* metric calculations
* biological validation
* report formatting

Those responsibilities belong to their respective packages.

---

# Workflow

A typical runner follows this sequence

```
Input annotation

↓

Parser

↓

Standard annotation object

↓

Metric calculation

↓

Metric results

↓

Report generation

↓

Output files
```

The runner is responsible for connecting these stages.

---

# Inputs

Runners typically receive

* annotation file paths
* command-line arguments
* configuration options
* output directories

---

# Outputs

A runner may produce

* report files
* log messages
* workflow summaries
* exit status

The runner itself should not perform calculations beyond simple orchestration.

---

# Design principles

A runner should

* remain small
* remain readable
* delegate work to specialised packages
* contain minimal business logic

A runner should answer the question

> "What steps should be executed?"

rather than

> "How is this metric calculated?"

---

# Example workflow

```text
Run Feature-wide QC

↓

Parse annotation

↓

Compute feature-wide metrics

↓

Generate JSON report

↓

Generate TSV report

↓

Exit
```

---

# Adding a new runner

When creating a new workflow

1. Choose the required parser.
2. Load annotation data.
3. Call one or more metric modules.
4. Pass results to reports.
5. Handle logging and configuration.
6. Add integration tests.

---

# Common mistakes

❌ Writing parsing logic inside the runner

❌ Implementing metrics directly in the runner

❌ Formatting reports manually

❌ Duplicating workflow code across runners

---

# Testing

Runner tests should verify

* correct workflow execution
* integration between packages
* expected outputs
* error handling

They should not duplicate parser or metric unit tests.

---

# Related packages

* `parsers/`
* `metrics/`
* `reports/`
* `config/`

The runner is the integration point for the Annotation QC framework. Its role is to coordinate the workflow while keeping each package independent and reusable.
