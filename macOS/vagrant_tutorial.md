# This is guidlines for working with [vagrant on macOS](https://developer.hashicorp.com/vagrant/docs/other/macos-catalina)
Main purpose is for the purpose of creating sif from docker containers

## Usage
* Enter directory with `Vagrantfile`
* Run `vagrant up` to startup vm
* Enter vm with `vagrant ssh`
* Kill vm with `vagrant halt`


## This is the Vagrantfile used to build env with singularity and docker
```
Vagrant.configure("2") do |config|
  config.vm.box = "singularityware/singularity-2.4"
  config.vm.hostname = "ubuntu"

  config.vm.provider :docker do |docker, override|
    override.vm.box = nil
    docker.image = "rofrano/vagrant-provider:ubuntu"
    docker.remains_running = true
    docker.has_ssh = true
    docker.privileged = true
    docker.volumes = ["/sys/fs/cgroup:/sys/fs/cgroup:rw"]
    docker.create_args = ["--cgroupns=host"]
    # Uncomment to force arm64 for testing images on Intel
    # docker.create_args = ["--platform=linux/arm64"]     
  end  
end
```
