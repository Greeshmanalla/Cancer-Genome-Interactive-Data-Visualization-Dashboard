import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

#Reading CSV file into dataframe
mutational_df = pd.read_csv('BRCA_Mutational.csv',low_memory=False)
clinical_df = pd.read_csv('BRCA_Clinical.csv')

# Define the initial list of irrelevant values
initial_irrelevant_values = [
    'synonymous_variant', 'incomplete_terminal_codon_variant',
    'intron_variant', 'mature miRNA_variant', 'inframe_deletion',
    'inframe_insertion', 'non_coding_transcript_exon_variant'
]

# Define additional irrelevant values
additional_irrelevant_values = [
    '3_prime_UTR_variant', '5_prime_UTR_variant',
    'coding_sequence_variant', 'downstream_gene_variant',
    'upstream_gene_variant'
]

# Combine both lists of irrelevant values
all_irrelevant_values = initial_irrelevant_values + additional_irrelevant_values

# Filter out rows where One_Consequence column contains irrelevant values
filtered_df = mutational_df[~mutational_df['One_Consequence'].isin(all_irrelevant_values)]

# Check for null values in the entire Clinical DataFrame
null_values_clinical = clinical_df.isnull().sum()

# Print the count of null values for each column
print("Null values in each column:")
print(null_values_clinical)


# Drop rows where all values are null in clinical df
clinical_df = clinical_df.dropna(how='all')

# Remove columns with entirely blank values in clinical df
cleaned_clinical_df = clinical_df.dropna(axis=1, how='all')


# Check for null values in the entire Mutational DataFrame
clinical_cleaned_null_values = cleaned_clinical_df.isnull().sum()

# Fill null values with 'NA'
filled_clinical_df = cleaned_clinical_df.fillna(value='NA')

#Just verifying null values still present or not?
# Check for null values
null_value = filled_clinical_df.isnull().sum()

# Check for null values in the entire Mutational DataFrame
null_values = filtered_df.isnull().sum()

# Print the count of null values for each column
print("Null values in each column:")
print(null_values)
# Drop rows where all values are null
filtered_df = filtered_df.dropna(how='all')


# Fill null values with 'NA'
filled_df = filtered_df.fillna(value='NA')

null_values_latest = filled_df.isnull().sum()
print(null_values_latest)

# Extract patient ID from Tumor_Sample_barcode and knowing unique no. of patients..

# Extract patient IDs
def extract_patient_id(barcode):
    parts = barcode.split('-')
    if len(parts) > 2:
        return parts[2]
    else:
        return None  # or some default value

filled_df['patient_id'] = filled_df['Tumor_Sample_Barcode'].apply(extract_patient_id)

# filled_df['patient_id'] = filled_df['Tumor_Sample_Barcode'].apply(lambda x: x.split('-')[2])

# Get the distinct number of rows count in the column 'Patient ID'
# Getting number of distinct patients..
unique_patients_count = filled_df['patient_id'].nunique()

print(f"Number of unique patients in mutational dataset: {unique_patients_count}")

# Drop un-necessary columns in clinical dataset..
columns_to_drop_clinical = ['days_to_death']
# Drop specified columns
filled_clinical_df.drop(columns=columns_to_drop_clinical, inplace=True)
# print(filled_clinical_df.columns)
#filled_clinical_df['Patient_Id'].nunique()

# Remove duplicates by keeping the first occurrence of both datasets
clinical_df = filled_clinical_df.drop_duplicates(subset='patient_id', keep='first')
# clinical_df = filled_clinical_df.drop_duplicates(subset='Patient_Id', keep='first')
brca_mutational_df = filled_df.drop_duplicates(subset='patient_id', keep='first')

# No. of unique patients in clinical dataset..
print(f"Number of unique patients in clinical dataset: {filled_clinical_df['patient_id'].nunique()}")
filled_clinical_df['patient_id'].isnull().sum()

# Merge the datasets on 'patient_id'
merged_df = pd.merge(clinical_df, brca_mutational_df, on='patient_id', how='outer')

##
#10:13AM with color blind safe colors and better readability..
import plotly.express as px
from dash import Dash, dcc, html
from dash.dependencies import Input, Output

# Load your datasets (assuming they are already preprocessed and cleaned)
df = merged_df
mutations_df = filled_df

# Initialize the Dash app
app = Dash(__name__)

app.layout = html.Div([
    html.H1("Cancer Genome Patient Data Dashboard"),

    html.Div([
        html.Label("Select Variable for Distribution:"),
        dcc.Dropdown(
            id='distribution-var',
            options=[{'label': col, 'value': col} for col in ['gender', 'race_list', 'ethnicity', 'tumor_tissue_site']],
            value='gender'
        )
    ]),

    dcc.Graph(id='distribution-plot'),

    html.Div([
        html.Label("Select X and Y variables for Scatter Plot:"),
        dcc.Dropdown(
            id='scatter-x-var',
            options=[{'label': col, 'value': col} for col in df.columns],
            value='days_to_last_followup'
        ),
        dcc.Dropdown(
            id='scatter-y-var',
            options=[{'label': col, 'value': col} for col in df.columns],
            value='age_at_initial_pathologic_diagnosis'
        )
    ]),

    dcc.Graph(id='scatter-plot'),

    html.H1("Gene Mutations Stacked Bar Chart"),
    dcc.Graph(id='stacked-bar-chart')
])

# Define a color-blind-safe palette
color_palette = px.colors.qualitative.Safe

@app.callback(
    Output('distribution-plot', 'figure'),
    Input('distribution-var', 'value')
)
def update_distribution_plot(selected_var):
    if selected_var not in df.columns:
        return {}
    fig = px.histogram(df, x=selected_var, title=f'Distribution of {selected_var}', color_discrete_sequence=color_palette)
    fig.update_layout(
        xaxis_title=selected_var.capitalize(),
        yaxis_title='Count',
        font=dict(size=14)
    )
    return fig

@app.callback(
    Output('scatter-plot', 'figure'),
    [Input('scatter-x-var', 'value'),
     Input('scatter-y-var', 'value')]
)
def update_scatter_plot(x_var, y_var):
    if x_var not in df.columns or y_var not in df.columns:
        return {}
    fig = px.scatter(df, x=x_var, y=y_var, title=f'{y_var} vs {x_var}', color_discrete_sequence=color_palette)
    fig.update_layout(
        xaxis_title=x_var.replace('_', ' ').capitalize(),
        yaxis_title=y_var.replace('_', ' ').capitalize(),
        font=dict(size=14)
    )
    return fig


@app.callback(
    Output('stacked-bar-chart', 'figure'),
    Input('stacked-bar-chart', 'id')
)
def update_stacked_bar_chart(_):
    mutation_counts = mutations_df.groupby(['Hugo_Symbol', 'One_Consequence']).size().reset_index(name='Mutation_Count')
    total_mutations = mutation_counts.groupby('Hugo_Symbol')['Mutation_Count'].sum().reset_index(name='Total_Mutations')
    mutation_counts = mutation_counts.merge(total_mutations, on='Hugo_Symbol')
    mutation_counts = mutation_counts.sort_values(by='Total_Mutations', ascending=False)
    gene_order = mutation_counts['Hugo_Symbol'].unique()

    fig = px.bar(
        mutation_counts,
        x='Hugo_Symbol',
        y='Mutation_Count',
        color='One_Consequence',
        title='Number of Mutations per Gene colored by Mutation Type',
        labels={'Mutation_Count': 'Number of Mutations'},
        category_orders={'Hugo_Symbol': gene_order},
        color_discrete_sequence=color_palette
    )
    fig.update_layout(
        xaxis_title="Gene",
        yaxis_title="Number of Mutations",
        legend_title="Mutation Type",
        font=dict(
            size=14
        )
    )
    return fig

if __name__ == '__main__':
    app.run_server(port=8060)
