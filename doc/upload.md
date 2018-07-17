# Data upload instructions

To upload the generated data to your instance of pretzel, you may use curl POST requests

## Retrieve your authentication token from pretzel front-end and record it 

```
TOKEN="your-token-string-here"
```

## Upload dataset (genome) definitions

```
time for F in *_genome.json; do 
  curl -X POST --header 'Content-Type: application/json' \
  --header 'Accept: application/json' -d @${F} \
  http://localhost:3000/api/Datasets/createComplete?access_token=${TOKEN} 
done
```

## Upload features (genes) definitions

```
time for F in *_annotation.json; do 
  curl -X POST --header 'Content-Type: application/json' \
  --header 'Accept: application/json' -d @${F} \
  http://localhost:3000/api/Datasets/createComplete?access_token=${TOKEN}
done
```


## Upload aliases 

```
for F in *_aliases.json; do 
  curl -X POST --header 'Content-Type: application/json' \
  --header 'Accept: application/json' -d @${F} \
  http://localhost:3000/api/Aliases/bulkCreate?access_token=${TOKEN}
done 
```