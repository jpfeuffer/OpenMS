void setBuildStatus(String message, String state) {
  step([
      $class: "GitHubCommitStatusSetter",
      reposSource: [$class: "ManuallyEnteredRepositorySource", url: "https://github.com/my-org/my-repo"],
      contextSource: [$class: "ManuallyEnteredCommitContextSource", context: "ci/jenkins/build-status"],
      errorHandlers: [[$class: "ChangingBuildStatusErrorHandler", result: "UNSTABLE"]],
      statusResultSource: [ $class: "ConditionalStatusResultSource", results: [[$class: "AnyBuildResult", message: message, state: state]] ]
  ]);
}

pipeline {
    agent none

    stages {
        stage("build and deploy on Windows and Linux") {
            parallel {
                stage("windows") {
                    agent {
                        label "win"
                    }
                    stages {
                        stage("build") {
                            steps {
                                sh '''#!bash
                                  echo "win" > win.txt
                                '''
                            }
                        }
                        stage("deploy") {
                            when {
                                branch "master"
                            }
                            steps {
                                sh '''#!bash
                                  cat win.txt
                                '''
                            }
                        }
                    }
                }

                stage("linux") {
                    agent {
                        label "openms&&ci-ready"
                    }
                    stages {
                        stage("configure") {
                            steps {
                                sh '''
                                mkdir -p bld && cd bld
                                cmake -DCMAKE_PREFIX_PATH='/usr/;/usr/local' -DOPENMS_CONTRIB_LIBS="/contrib-build/" -DBOOST_USE_STATIC=OFF -DHAS_XSERVER=Off ../OpenMS
                                '''
                            }
                        }
                        stage("build") {
                            steps {
                                sh '''
                                cd bld
                                cmake --build . --target OpenMS --config Debug -j6
                                '''
                            }
                        }
                        stage("test") {
                            steps {
                                sh '''
                                cd bld
                                ctest -D Continuous -T test --group "pyopenms"
                                '''
                            }
                        }
                        stage("deploy") {
                             when {
                                 branch "nightly"
                             }
                             steps {
                                sh 'echo "I will deploy"'
                            }
                        }
                    }
                }
            }
        }
    }
}
