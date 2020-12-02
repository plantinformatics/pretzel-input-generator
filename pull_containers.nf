#!/usr/bin/env nextflow

import nextflow.util.Escape
import nextflow.container.SingularityCache

def containers = []
session.getConfig().process.each {k, v ->
  if((k.startsWith('withLabel:') || k.startsWith('withName:')) && v.containsKey('container')) {
    println "$k -> $v.container"
    containers << v.container
  }
}

SingularityCache scache = new SingularityCache() //to get NF-consitent image file names

process pull_container {
  container null
  tag { remote }
  maxForks 1
  storeDir "${params.singularitydir}"
  echo true

input:
  val(remote) from Channel.from(containers)

output:
  file(img)

script:
img = scache.simpleName(remote)
"""
singularity pull --name ${img} docker://${Escape.path(remote)}
"""
}