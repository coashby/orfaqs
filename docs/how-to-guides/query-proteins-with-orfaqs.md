# Load, Query, and Export Proteins with the ORFaqs CLI
The ORFaqs Protein Query tool (or Query tool for short) allows you to load and interact with protein data using query using SQL syntax. It also supports exporting stored protein data and query results.

The CLI allows you to:
- Load proteins data into a local relational database
- Query stored protein data
- Export stored protein data and queries
- Manage workspaces and tables

**Prerequisites**

If you are using the CLI from source, make sure you complete the items below before proceeding with this tutorial.
- [Setting Up Your Dev Environment](./../quick-start-quides/setting-up-your-dev-environment.md)
----


## Load Proteins into a Local Database

The `load-proteins` command helps you load protein data found using the ORFaqs
Protein Discovery tool, or reference proteins from FASTA proteome files. The
tool will automatically detect the type of data being loaded. By default, any
protein data loaded is placed in a table called `proteins` within the workspace
`orfaqs_protein_query`.

**Basic Usage**
```bash
orfaqs_protein_query load-proteins PROTEINS [OPTIONS]
```

When loading **discovered proteins**:
- Specify an individual `CSV` results file, or a directory of result files.
- If the path is a directory path, all discovered protein files within that path
are loaded (all other files are ignored).

When loading **reference proteins**:
- Specify an individual `FASTA` protein file, or a directory of `FASTA` protein
files.
- If the path is a directory path, all `FASTA` files within that path are
loaded.

### Set Workspace and Table Names
Use the `--workspace` and `--table-name` options to change the defaults used in
assigning the workspace a table name respectively.
```bash
orfaqs_protein_query load-proteins PROTEINS --workspace WORKSPACE --table TABLE
```

## Run Queries

The `query` command allows you to run SQL `SELECT` queries against proteins
stored in your workspace database. You can filter, aggregate, and analyze
protein data using standard SQL syntax.

**Basic Usage**
```bash
orfaqs_protein_query query WORKSPACE [QUERY] [OPTIONS]
```

Query results are then returned to the current terminal using an interactive
display.

### Display the Number of Proteins within a Table
```bash
orfaqs_protein_query query WORKSPACE "SELECT COUNT(*) FROM table"
```

### Display All Proteins Matching a Length Criteria Sorted by Length
```bash
orfaqs_protein_query query WORKSPACE "SELECT * FROM table WHERE protein_length >= 10 AND protein_length <= 50 ORDER BY protein_length ASC"
```

### Display the Number of Reference Proteins that Overlap with the Set of Discovered Proteins
```bash
orfaqs_protein_query query WORKSPACE "SELECT COUNT(rp.*) FROM reference_proteins_table rp WHERE EXISTS (SELECT 1 FROM discovered_proteins_table dp WHERE rp.protein = dp.protein)"
```

### Export Query Results
The `--export` option exports the results of the query to the default directory
using the default export format (*CSV*). You can set the export file path and
format by using the `--export-path` and `--export-format` options respectively.

#### Export a Query Using the Default Options
```bash
orfaqs_protein_query query WORKSPACE [QUERY] --export
```

#### Export a Query as an Excel File
```bash
orfaqs_protein_query query WORKSPACE [QUERY] --export-format xlsx
```

#### Export a Query in FASTA format
```bash
orfaqs_protein_query query WORKSPACE [QUERY] --export-format fasta
```