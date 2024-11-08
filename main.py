import numpy as np
import pandas as pd


#Function that manage the duplations in the Array
def remove_duplicates_2d_array(array):
    seen = set()
    unique_array = []
    for row in array:
        row_tuple = tuple(row)
        if row_tuple not in seen:
            seen.add(row_tuple)
            unique_array.append(row)
    return unique_array
#Function that sort based on Probabilities 
def sort_2d_array_by_second_element(array):
    return sorted(array, key=lambda x: x[1], reverse=True)

#Highly abundant proteins's IDs
prot= [
    "P02768", "P02787", "P00738", "P01876-1", "P02671", "P02675", "P01024", 
    "P02679", "P02647", "P02790", "P02652", "P02774", "P01023", "P01860", 
    "P02765", "P08603", "P01871-1", "P00450", "P00751", "P01857-1", 
    "P01042-2", "P0DOX7", "P02763", "P00747", "P04217"
]



tab = list()

#Loading Data
data = pd.read_csv('Supp_Data.csv')

#Finding Min and Max vallues of the Log2FC
f = data['Log2FC']
mini = round(min(f))*10
maxi = round(max(f))*10

# z and j are the theresholds 
for z in range(1,maxi,5):
    for j in range(mini,1,5):
        temp_i = z/10
        temp_j = j/10
        #defining effect based on the new theresholds
        data['result'] = data.apply(lambda row: 1 if ((row['Accession']) in prot and row['Log2FC']) <= temp_j or ((row['Accession']) not in prot and row['Log2FC']) >= temp_i else 0, axis=1)
        df2 = data.loc[:, ['description', 'Group','concentration','result']]
        df = df2.copy()

        col = df.columns
        cause = []

        #Causal Analysis
        for index, row in df.iterrows():
            # Calculating the actual world
            factual = df[(df[col[0]] == row[0]) & (df[col[1]] == row[1]) & (df[col[2]] == row[2]) & (df[col[3]] == 1)]
            # Calculating the counterfactual world
            counteractual = df[(df[col[0]] == row[0]) & (df[col[1]] != row[1]) & (df[col[2]] == row[2]) &(df[col[3]] == 1)]
            actual = df[(df[col[1]] == row[1]) & (df[col[3]] == 1)]
            all_times = df[df[col[1]] == row[1]]
            #checking whether probability of happening the actual world is greater than counterfactual world (counterfactual reasoning)
            if len(factual)> len(counteractual):
                #adding the possible cause and probability of it
                cause.append([row[1],len(actual)/len(all_times)])
        #removinf duplicates
        uni = remove_duplicates_2d_array(cause)
        #sort
        uni_sorted = sort_2d_array_by_second_element(uni)

        #Extarct top Three highest probability
        if len(uni_sorted)>= 3 :
            l = uni_sorted[:3]
            tab.append([temp_i,temp_j,l[0],l[1],l[2]])
        elif len(uni_sorted)== 2:
            tab.append([temp_i,temp_j,uni_sorted[0],uni_sorted[1],None])
        elif len(uni_sorted)== 1:
            tab.append([temp_i,temp_j,uni_sorted[0],None,None])
        else:
            break

fin_tab = list()
for i in tab:
    st = lsi()
    for z in range(2,5):
        if i[z]==None:
            st.append('None')
        else:
            s = i[z][0] +' with prob '+ str(i[z][1])
            st.append(s)

    fin_tab.append([i[0],i[1],st[0],st[1],st[2]])

df = pd.DataFrame(fin_tab, columns=['high threshold', 'low threshold', 'cause 1', 'cause 2', 'cause 3'])

df.to_csv('result.csv',index=False)


