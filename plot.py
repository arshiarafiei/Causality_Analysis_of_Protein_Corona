import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('result.csv')


df['result'] = df['casue 1'].apply(lambda x: x.split(' with')[0])

df1 = df.loc[:, ['high threshold', 'low threshold','result']]
df = df1.copy()


df['result_code'] = df['result'].astype('category').cat.codes


pivot_table = df.pivot(index="low threshold", columns="high threshold", values="result_code")


plt.figure(figsize=(30, 30))
sns.heatmap(pivot_table, annot=df.pivot(index="low threshold", columns="high threshold", values="result"), fmt='', cmap='viridis')
plt.title('Heatmap of causes')
plt.savefig('heatmap.png')
plt.show()