import numpy as np
import os
import pdb
import json
import glob
from nltk.corpus   import stopwords
from nltk.tokenize import word_tokenize
from nltk          import WordNetLemmatizer

np.random.seed(42)

PATH      = "./dataset/bbc/2020/"
MONTHS    = [str(i) for i in range(1,6)]
STOPWORDS = set(stopwords.words('english'))

def argsort(seq):
    return [i for (v, i) in sorted((v, i) for (i, v) in enumerate(seq))]

def main():

    ## construct dataset

    lemma = WordNetLemmatizer()
    content_per_month = {}
    
    word2idx = {}
    idx2word = []
    corpus   = []

    for month in MONTHS:
        files      = glob.glob(os.path.join(PATH, month) + '/**/articles', recursive=True)
        days       = [int(f.split('/')[5]) for f in files]
        sorted_idx = argsort(days)

        content_per_month[int(month)] = []

        for idx in sorted_idx:
            with open(files[idx], 'r') as f:
                articles = json.load(f)
                f.close()

            data   = articles['articles']
            sample = np.random.choice(len(data), 10, replace=False)
            
            content_per_day = []
            title_per_day   = []

            n = 1
            for sample_idx in sample:
                
                if n > 5:
                    break

                r_content = data[sample_idx]['content']            
                r_title   = data[sample_idx]['title']
                
                content = word_tokenize(r_content)
                title   = word_tokenize(r_title)
                
                content = [w.lower() for w in content] 
                title   = [w.lower() for w in title] 
                
                content = [w for w in content if not w in STOPWORDS]
                title   = [w for w in title if not w in STOPWORDS] 
                
                content = [w for w in content if w.isalpha()] 
                title   = [w for w in title if w.isalpha()] 
                
                content = [lemma.lemmatize(w, pos="v") for w in content] 
                title   = [lemma.lemmatize(w, pos="v") for w in title] 
                
                content = [lemma.lemmatize(w, pos="n") for w in content] 
                title   = [lemma.lemmatize(w, pos="n") for w in title] 

                if len(content) < 100:
                    continue
                else:
                    n += 1
                
                content = np.random.choice(content, 100, replace=False).tolist() # sample
                                
                content_per_day += content
                title_per_day   += title

                corpus += content

            content_per_month[int(month)].append(content_per_day)
    
    ## construct vocabulary

    n = 0
    for token in corpus:
        if token not in word2idx.keys():
            word2idx[token] = n
            idx2word.append(token)
            n += 1
    
    ## token to index

    content_per_month_idx = {}

    for month in range(1, 6):
        content_per_month_idx[month] = []
        for article_day in content_per_month[month]:
            article_day_idx = []
            for token in article_day:
                article_day_idx.append(word2idx[token])                               
            content_per_month_idx[month].append(article_day_idx)

    ## dump files

    dump = []
    for month in range(1, 6):
        for article in content_per_month_idx[month]:
            dump.append(np.array(article))
    
    dump     = np.array(dump)
    idx2word = np.array(idx2word)

    np.save('dump.npy', dump)
    np.save('idx2word.npy', idx2word)

    print("Fin.")

    
        
                
if __name__=="__main__":
    main()


