Vagrant.configure('2') do |config|
  config.vm.box = 'ubuntu'

  config.vm.provider :aws do |aws, override|
    aws.region = "us-east-1"
    aws.security_groups = "quick-start-1" 
    
    # 8 cores @ $0.82/h
    aws.instance_type = "m3.2xlarge"
    # 32 cores @ $2.40/h
    # aws.instance_type = "cc2.8xlarge"

    aws.tags = {
      'Name' => 'BMRF-Compute-Node',
    }
   
    override.ssh.username = "ubuntu"
  end

  config.vm.provision :shell,
     :path => "devops/provision.sh"

end
