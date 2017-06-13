echo "##Author: Sean Liu" > README.md
date | xargs -I{} echo "##TIME:" {} >> README.md
echo "
我和你不再联系
希望你不要介意
要怪就怪当初没在一起
而你对现在也比较满意
所以我留下来也没有道理
我和你断了联系
不代表我不想你
走到哪里还是会有惦记
而我也开始试着去忘记
抹去我们过去的放弃的所有交集
" >> README.md

git add `git ls-files -o`
git add `git ls-files -m`
git rm `git ls-files -d`
date | xargs -I{} git commit -m {}
git push origin master
