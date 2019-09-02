# Data upload instructions

To upload the generated data to your instance of pretzel, you may use curl POST requests. Make sure these are executed _in order_.

## 1. Retrieve your authentication token from pretzel front-end and record it.

In firefox, open menu and select *Web Developer*, *Network* then refresh, select one of the requests, and click  *Cookies* to find your token string.

```
TOKEN="your-token-string-here"
```

Point to your instance of pretzel, for example:

```
SRV="http://localhost:3000"
```

## 2. Upload dataset (genome) definitions

```
for F in *_genome.json; do
  curl -X POST --header 'Content-Type: application/json' \
  --header 'Accept: application/json' -d @${F} \
  "${SRV}/api/Datasets/createComplete?access_token=${TOKEN}"
  echo
done
```

## 3. Upload marker sequences and features (genes) definitions (compressed)

```
for F in *_{markers,annotation}.json.gz; do
  curl -X POST --header 'Content-Type: application/json' \
  --header 'Accept: application/json' -H'Content-Encoding: gzip' \
  --data-binary @${F} \
  "${SRV}/api/Datasets/createComplete?access_token=${TOKEN}"
  echo
done
```

## 4. Upload aliases (compressed)

```
for F in *_aliases.json.gz; do
  echo -ne "\n${F}\t"
  curl -X POST --header 'Content-Type: application/json'   \
  --header 'Accept: application/json' -H'Content-Encoding: gzip' \
  --data-binary @${F}   \
  "${SRV}/api/Aliases/bulkCreate?access_token=${TOKEN}"
done
```

## Troubleshooting

If there are too many aliases for your isntance of pretzel to handle, leading to out of memory errors, there are several things which may help:

* reduce the number of aliases by increasing [filtering stringency in your conf/input.config](https://github.com/plantinformatics/pretzel-input-generator/blob/v1.0/conf/input.config#L4-L9)
* split your alias file and upload in chunks
* re-run `node` with more memory, e.g. (`--max-old-space-size=8192`).


# Uploading uncompressed data

## Upload features (genes) definitions

```
for F in *_annotation.json; do
  curl -X POST --header 'Content-Type: application/json' \
  --header 'Accept: application/json' -d @${F} \
  "${SRV}/api/Datasets/createComplete?access_token=${TOKEN}"
done
```


## Upload aliases (plain text)

```
for F in *_aliases.json; do
  curl -X POST --header 'Content-Type: application/json' \
  --header 'Accept: application/json' -d @${F} \
  "${SRV}/api/Aliases/bulkCreate?access_token=${TOKEN}"
done
```

## Deleting existing data

Before re-uploading of updated datasets, delete the existing ones (requires [keyFromJSON.py](https://github.com/plantinformatics/pretzel-input-generator/blob/master/bin/keyFromJSON.py) for extracting dataset names):

```
for name in $(./keyFromJSON.py JSON/*_{annotation,genome}.json*); do
  echo -ne  "\n\nTrying to delete: ${name}"
  curl -X DELETE --header 'Accept: application/json' \
  "${SRV}/api/Datasets/${name}?access_token=${TOKEN}"
done; echo
```


