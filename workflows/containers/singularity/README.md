# Singularity containers for the nextflow workflow

This document described how singularity containers are build from .def files (described 
in more detail [here](https://github.com/bhattlab/bhattlab_workflows/tree/master/containers/preprocessing))
and how those are then published to the Sylab could, to be available independent of the
file system.

## Preparations

You will need an account with the [Sylabs could](https://cloud.sylabs.io/). You can just log in with your
GitHub credentials and generate a new account.

#### Generate an access token

In the [dashboard](https://cloud.sylabs.io/dashboard), you can generate an access token. Save this somewhere
because you will not be able to see it again.

#### Configuration file

Next, you can configure your HPC system (in our case, the SCG cluster) to use the remote capabilities
of the Sylabs cloud. To do so, generate a configuration file:
```
touch ~/.singularity/remote.yaml
```
and add these information:
```
Active: SylabsCloud
Remotes:
  SylabsCloud:
    URI: cloud.sylabs.io
    Token: <your access token>
    System: true
    Exclusive: false
```

#### Generate a keypair

You will also need a keypair to sign your container. You can do this via:
```
$ singularity key newpair
Enter your name (e.g., John Doe) : John Doe
Enter your email address (e.g., john.doe@example.com) : John.Doe@example.com
Enter optional comment (e.g., development keys) : scs demo keys
Enter a passphrase :
Would you like to push it to the keystore? [Y,n] y
Generating Entity and OpenPGP Key Pair... done
Key successfully pushed to keystore
```

## Build a container

Now you can use the remote build capacities of the Syslabs cloud to build the container that you have
previously generated via Docker 
(see [these instructions](https://github.com/bhattlab/bhattlab_workflows/tree/master/containers/preprocessing)).

```
singularity build --remote <image> <def>
```
For the preprocessing container, it would look something like that:
```
singularity build --remote bhattlab-micromamba-preprocessing.img preprocessing.def
```

## Sign a container

To publish the container, it will have to be signed. You can sign it with this command:
```
$ singularity sign <image>
Enter encryption passphrase :
```
This container can then also be verified via singularity:
```
singularity verify <image>
```

## Publish via Sylabs

Finally, you can upload your container to the Sylabs cloud, where it then can be published so that it is
available for other people as well:

```
singularity push <image> library://<user>/<repo>/<image>:<tag>
```


For me (wirbel), and for the preprocessing container, the command would look like this:
```
singularity push bhattlab-micromamba-preprocessing.img library://wirbel/bhattlab/bhattlab-micromamba-preprocessing.sif:v0.2
```


