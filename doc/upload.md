# Data upload instructions

To upload the generated data to your instance of pretzel, you may use curl POST requests

## Retrieve your authentication token from pretzel front-end and record it.

In firefox, open menu and select *Web Developer*, *Network* then refresh, select on of the requests, and click  *Cookies*

```
TOKEN="your-token-string-here"
```

Point to your instance of pretzel

```
SRV="https://localhost:3000"
```

## Upload dataset (genome) definitions

```
for F in *_genome.json; do
  curl -X POST --header 'Content-Type: application/json' \
  --header 'Accept: application/json' -d @${F} \
  "${SRV}/api/Datasets/createComplete?access_token=${TOKEN}"
done
```

## Upload features (genes) definitions

```
for F in *_annotation.json; do
  curl -X POST --header 'Content-Type: application/json' \
  --header 'Accept: application/json' -d @${F} \
  "${SRV}/api/Datasets/createComplete?access_token=${TOKEN}"
done
```

## Upload aliases (compressed)

```
for F in *_aliases.json.gz; do
  curl -X POST --header 'Content-Type: application/json'   \
  --header 'Accept: application/json' -H'Content-Encoding: gzip' \
  --data-binary @${F}   \
  "${SRV}/api/Aliases/bulkCreate?access_token=${TOKEN}"
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

## Troubleshooting

If there are too many aliases for your isntance of pretzel to handle, leading to out of memory errors, there are several things which may help:

* reduce the number of aliases by incresing [filtering stringency in your conf/input.config](https://github.com/plantinformatics/pretzel-input-generator/blob/d4e7c88776c5f9c4ab6f9d50adcf49bd36cf6f81/conf/input.config#L4-L9)
* split your alias file and upload in chunks
* re-run `node` with more memory, e.g. (`--max-old-space-size=8192`).
