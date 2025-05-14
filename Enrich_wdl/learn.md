
## 使用eggnog-mapper拿到建背景库准备文件
[最好的建库教程！模式植物构建orgDb数据库 | 以org.Slycompersicum.eg.db为例](https://mp.weixin.qq.com/s/b8OrDKJJGdXwF9B1C7l6zg)

1. 该物种的蛋白质序列，最好蛋白序列的名称和单细胞测序数据的基因名是一致的
2. 提交到[eggnog-mapper](http://eggnog-mapper.embl.de/)
选中蛋白质序列→填写邮箱→提交后点击邮箱信息→开始job→下载xlsx文件

Check部分是检查输入基因集是否在建好的背景库中，很大概率不是完全一致，需要个性化的添加后缀等，另外下一个`task: enrich`要求输入的csv必须要这三列`gene_id`, `cluster`, `p_val_adj`，所以在`check`部分可以编辑代码使其满足该结构