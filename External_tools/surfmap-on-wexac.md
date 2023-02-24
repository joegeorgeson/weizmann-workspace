

For pulling the docker and initial setup;
```
ssh docker1
docker pull nchenche/surfmap:v1.5
docker images |grep nchenche/surfmap:v1.5
docker tag 58e7e4fd08d5 ops:5000/surfmapv15 #the tag here should be unique
docker push ops:5000/surfmapv15
singularity pull docker://nchenche/surfmap:v1.5
```

More to come
