import pandas as pd
import matplotlib.pyplot as plt

# Creating a DataFrame from the provided data
data = {
    'Day': ['Day1', 'Day1', 'Day2', 'Day2', 'Day3', 'Day3', 'Day4', 'Day4', 'Day5', 'Day5', 'Day6', 'Day6', 'Day7', 'Day7'],
    'Model': ['HR3a', 'HR1', 'HR3a', 'HR1', 'HR3a', 'HR1', 'HR3a', 'HR1', 'HR3a', 'HR1', 'HR3a', 'HR1', 'HR3a', 'HR1'],
    'Bias': [0.0075, -0.1248, 0.0214, -0.1118, 0.0339, -0.1017, 0.0442, -0.0950, 0.05098, -0.0908, 0.0616, -0.0822, 0.06811, -0.0825],
    'RMSE': [0.2830, 0.1843, 0.2803, 0.1614, 0.2843, 0.1596, 0.2938, 0.1747, 0.3163, 0.2002, 0.3372, 0.2368, 0.3636, 0.2793]
}

df = pd.DataFrame(data)

# Plotting the bias for both models over the days
plt.figure(figsize=(12, 10))

# Plotting Bias
plt.subplot(2, 1, 1)  # 2 rows, 1 column, 1st subplot
colors = {'HR3a': 'blue', 'HR1': 'red'}  # Dictionary to hold the colors for each model
for model in df['Model'].unique():
    model_data = df[df['Model'] == model]
    plt.plot(model_data['Day'], model_data['Bias'], marker='o', color=colors[model], label=f'{model} Bias')

plt.xlabel('Days')
plt.ylabel('Bias')
plt.title('Bias of Models HR3a and HR1 Over Days')
plt.legend()
plt.grid(True)
plt.xticks(rotation=45)

# Plotting RMSE
plt.subplot(2, 1, 2)  # 2 rows, 1 column, 2nd subplot
for model in df['Model'].unique():
    model_data = df[df['Model'] == model]
    plt.plot(model_data['Day'], model_data['RMSE'], marker='o', color=colors[model], label=f'{model} RMSE')

plt.xlabel('Days')
plt.ylabel('RMSE')
plt.title('RMSE of Models HR3a and HR1 Over Days')
plt.legend()
plt.grid(True)
plt.xticks(rotation=45)

plt.tight_layout()
#plt.show()

# Saving the plot as a PNG file
plt.savefig('./Bias_RMSE_plot_summer.png')  # Adjusted path for saving in this environment

