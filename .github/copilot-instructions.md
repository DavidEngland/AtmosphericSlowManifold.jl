# Figure and Table Generation Guidelines

When requested to generate, convert, or format tables or figures for reports:

1. Always provide dual-format outputs (Markdown + LaTeX) unless explicitly asked for only one.

2. Table Standards:
   - Markdown: Use GitHub-Flavored Markdown (GFM) with explicit column alignment (`| :--- | :---: | ---: |`).
   - LaTeX: Use the `table` environment paired with `booktabs` (`\toprule`, `\midrule`, `\bottomrule`). Never use vertical lines (`|`) in professional LaTeX tables. Include `\caption{}` and `\label{tbl:...}`.

3. Figure Standards:
   - Markdown: Use standard image syntax `![Caption](path/to/image.png)` or embedded Mermaid.js diagrams (```mermaid ... ```) for schematics.
   - LaTeX: Use standard `\begin{figure}[htbp] \centering ... \end{figure}` wrappers with `\includegraphics[width=\linewidth]{...}`, `\caption{}`, and `\label{fig:...}`. For pure programmatic diagrams, provide TikZ code inside `\begin{tikzpicture}`.

4. Label Conventions:
   - Tables: `tbl:short-description` or `tab:short-description`
   - Figures: `fig:short-description`

## Interactive Prompts for VS Code Chat

Use these prompt templates when generating report artifacts.

### Prompt 1: Generating Dual-Format Tables from Raw Data or Descriptions

```text
Generate both a Markdown table and a publication-ready LaTeX table from the following data/description.

Data / Description:
[Insert raw CSV, JSON, bullet points, or text description here]

Requirements:
1. Markdown: Format with GFM syntax and clean column alignment. Add a bold title above the table.
2. LaTeX:
   - Wrap in a floating `table` environment with `[htbp]`.
   - Use `booktabs` package syntax (`\toprule`, `\midrule`, `\bottomrule`).
   - Center numbers and align text left.
   - Include a descriptive `\caption{}` and `\label{tab:...}`.
3. If applicable, add a footnote for abbreviations or mathematical symbols.
```

### Prompt 2: Generating Dual-Format Figures & Diagrams

```text
Create figure references for both Markdown and LaTeX for the following visualization/diagram.

Figure Details:
- Title/Concept: [e.g., Atmospheric Boundary Layer Profiling Architecture]
- File Path / Asset: [e.g., media/diagrams/abl_profile.png]
- Caption: [Insert descriptive caption]
- Label: [e.g., fig:abl-profile]

Requirements:
1. Markdown Output:
   - Standard image tag with caption.
   - Equivalent Mermaid.js diagram code block if this is a flowchart, state machine, or block architecture.
2. LaTeX Output:
   - Floating `figure` environment (`[htbp]`).
   - Centered alignment with `\linewidth` scaling.
   - Descriptive `\caption{}` and `\label{fig:...}`.
   - If applicable, provide standalone TikZ/PGFPlots code for native vector compilation.
```

### Prompt 3: Format Conversion (Markdown <-> LaTeX)

```text
Convert the following [Markdown / LaTeX] [table / figure] into the opposite format ([LaTeX / Markdown]) while maintaining structural fidelity.

Input Code:
[Paste snippet here]

Requirements:
- Preserve all captions, label identifiers, alignments, headers, and footnotes.
- For LaTeX output: ensure `booktabs` guidelines are enforced (no vertical lines).
- For Markdown output: convert LaTeX cross-references (`\ref{...}`) to clean inline references.
```

## Reference Output Templates

### 1) Dual Table Format Reference

#### Markdown Target

**Table 1: Physical Parameters for Boundary Layer Discretization**

| Parameter | Symbol | Range / Default | Units |
| :--- | :---: | ---: | :--- |
| Reference Boundary Height | $z_0$ | $0.0$ | $\mathrm{m}$ |
| Domain Top Boundary | $H$ | $1000.0$ | $\mathrm{m}$ |
| Gegenbauer Modes | $N$ | $6$ | --- |
| Viscous Scale | $\nu$ | $10^{-5}$ | $\mathrm{m^2 \, s^{-1}}$ |

#### LaTeX Target

```latex
\begin{table}[htbp]
  \centering
  \caption{Physical Parameters for Boundary Layer Discretization}
  \label{tab:physical-parameters}
  \begin{tabular}{lccc}
    \toprule
    \textbf{Parameter} & \textbf{Symbol} & \textbf{Range / Default} & \textbf{Units} \\
    \midrule
    Reference Boundary Height & $z_0$ & $0.0$ & $\mathrm{m}$ \\
    Domain Top Boundary & $H$ & $1000.0$ & $\mathrm{m}$ \\
    Gegenbauer Modes & $N$ & $6$ & --- \\
    Viscous Scale & $\nu$ & $10^{-5}$ & $\mathrm{m^2 \, s^{-1}}$ \\
    \bottomrule
  \end{tabular}
\end{table}
```

### 2) Dual Figure Format Reference

#### Markdown Target

*Figure 1: Atmospheric Slow Manifold Execution Pipeline.*

```markdown
![Atmospheric Slow Manifold Execution Pipeline](figures/pipeline_architecture.png)
```

#### LaTeX Target

```latex
\begin{figure}[htbp]
  \centering
  \includegraphics[width=0.85\linewidth]{figures/pipeline_architecture.png}
  \caption{Atmospheric Slow Manifold Execution Pipeline.}
  \label{fig:pipeline-architecture}
\end{figure}
```
