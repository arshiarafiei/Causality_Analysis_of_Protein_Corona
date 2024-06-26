import numpy as np
import pandas as pd



def remove_duplicates_2d_array(array):
    
    seen = set()
    unique_array = []
    for row in array:
        row_tuple = tuple(row)
        if row_tuple not in seen:
            seen.add(row_tuple)
            unique_array.append(row)
    return unique_array

def sort_2d_array_by_second_element(array):
    return sorted(array, key=lambda x: x[1], reverse=True)


prot = [
    ("P02768", "Albumin [OS=Homo sapiens]"),
    ("P0DOX5", "Immunoglobulin gamma-1 heavy chain [OS=Homo sapiens]"),
    ("P0DOX2", "Immunoglobulin alpha-2 heavy chain [OS=Homo sapiens]"),
    ("P0DOX6", "Immunoglobulin mu heavy chain [OS=Homo sapiens]"),
    ("P02675", "Fibrinogen beta chain [OS=Homo sapiens]"),
    ("P01009", "Alpha-1-antitrypsin [OS=Homo sapiens]"),
    ("P01023", "Alpha-2-macroglobulin [OS=Homo sapiens]"),
    ("P00738", "Haptoglobin [OS=Homo sapiens]"),
    ("P02786", "Transferrin receptor protein 1 [OS=Homo sapiens]"),
    ("P02647", "Apolipoprotein A-I [OS=Homo sapiens]"),
    ("P02652", "Apolipoprotein A-II [OS=Homo sapiens]"),
    ("P04114", "Apolipoprotein B-100 [OS=Homo sapiens]"),
    ("P01024", "Complement C3 [OS=Homo sapiens]"),
    ("P0C0L4", "Complement C4-A [OS=Homo sapiens]"),
    ("P00450", "Ceruloplasmin [OS=Homo sapiens]"),
    ("P02766", "Transthyretin [OS=Homo sapiens]"),
    ("P02763", "Alpha-1-acid glycoprotein 1 [OS=Homo sapiens]"),
    ("P00747", "Plasminogen [OS=Homo sapiens]"),
    ("P01011", "Alpha-1-antichymotrypsin [OS=Homo sapiens]"),
    ("P06396", "Gelsolin [OS=Homo sapiens]"),
    ("P69905", "Hemoglobin subunit alpha [OS=Homo sapiens]"),
    ("P08519", "Apolipoprotein(a) [OS=Homo sapiens]"),
    ("P02679", "Fibrinogen gamma chain [OS=Homo sapiens]")
]

tab = []

data = pd.read_csv('Supp_Data.csv')

f = data['Log2FC']
mini = round(min(f))*10
maxi = round(max(f))*10

for z in range(1,maxi,5):
    for j in range(mini,1,5):
        temp_i = z/10
        temp_j = j/10
        data['result'] = data.apply(lambda row: 1 if ((row['Accession'], row['description']) in prot and row['Log2FC']) <= temp_j or ((row['Accession'], row['description']) not in prot and row['Log2FC']) >= temp_i else 0, axis=1)
        df2 = data.loc[:, ['description', 'Group','concentration','result']]
        df = df2.copy()

        df2 = df[df['result'] == 1]
        col = df.columns
        cause = []
        for index, row in df2.iterrows():
            for i in range(len(row)-1):
                if i == 1:
                    df3 = df[(df[col[0]] == row[0]) & (df[col[1]] != row[1]) & (df[col[2]] == row[2]) & (df[col[3]] == 1)]
                    df4 = df[(df[col[0]] == row[0]) & (df[col[1]] != row[1]) & (df[col[2]] == row[2]) &(df[col[3]] == 0)]
                    counteractual = df[(df[col[1]] == row[1]) & (df[col[3]] == 0)]
                    factual = df[(df[col[1]] == row[1]) & (df[col[3]] == 1)]
                    all = df[df[col[1]] == row[1]]

                    if len(df3)+len(df4) == 0:
                        continue
                    else:
                        if len(df3)/(len(df3)+len(df4)) ==1:
                            # print(row[1])
                            # print(actual)
                            continue
                        else:
                            cause.append([row[i],len(factual)/len(all)])
                else:
                    continue
        uni = remove_duplicates_2d_array(cause)

        uni_sorted = sort_2d_array_by_second_element(uni)

        if len(uni_sorted)>= 3 :
            l = uni_sorted[:3]

            tab.append([temp_i,temp_j,l[0],l[1],l[2]])
        elif len(uni_sorted)== 2:
            tab.append([temp_i,temp_j,uni_sorted[0],uni_sorted[1],None])
        elif len(uni_sorted)== 1:
            tab.append([temp_i,temp_j,uni_sorted[0],None,None])
        else:
            break

fin_tab = []
for i in tab:
    st = []
    for z in range(2,5):
        if i[z]==None:
            st.append('None')
        else:
            s = i[z][0] +' with prob '+ str(i[z][1])
            st.append(s)

    fin_tab.append([i[0],i[1],st[0],st[1],st[2]])

df = pd.DataFrame(fin_tab, columns=['high threshold', 'low threshold', 'casue 1', 'cause 2', 'cause 3'])

df.to_csv('result.csv',index=False)
