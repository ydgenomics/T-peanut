Description: 拿到eggnog-mapper的xlsx文件后来建库，建库之后就可以跑富集分析，在跑之前之前保证库和目标基因集的名称是有较高覆盖度的！
Main task in WDL: 
`check`部分是检查输入基因集是否在建好的背景库中，很大概率不是完全一致，需要个性化的添加后缀等，另外下一个`task: enrich`要求输入的csv必须要这三列`gene_id`, `cluster`, `p_val_adj`，所以在`check`部分可以编辑代码使其满足该结构

## 使用eggnog-mapper拿到建背景库准备文件
[最好的建库教程！模式植物构建orgDb数据库 | 以org.Slycompersicum.eg.db为例](https://mp.weixin.qq.com/s/b8OrDKJJGdXwF9B1C7l6zg)

1. 该物种的蛋白质序列，最好蛋白序列的名称和单细胞测序数据的基因名是一致的  
    <img src="https://cloud.stomics.tech/workflow/api/static/68243e453bc6846f8fd1a1bc.png" alt="Description" width="500"/>

2. 提交到[eggnog-mapper](http://eggnog-mapper.embl.de/)  
    选中蛋白质序列→填写邮箱→提交submit→点击邮箱信息`Click to manage your job`→开始`Start job`→下载xlsx文件  
    <img src="https://cloud.stomics.tech/workflow/api/static/68243e9a06448a03866b11fe.png" alt="Description" width="500"/>
    <img src="https://cloud.stomics.tech/workflow/api/static/68243efb1f1f79cfb487e445.png" alt="Description" width="500"/>
    <img src="https://cloud.stomics.tech/workflow/api/static/68243f5d2b30ea192f5b886d.png" alt="Description" width="500"/>
    <img src="https://cloud.stomics.tech/workflow/api/static/682444d22b30ea192f5b8875.png" alt="Description" width="500"/>
    <img src="https://cloud.stomics.tech/workflow/api/static/682444f03bc6846f8fd1a1c1.png" alt="Description" width="500"/>