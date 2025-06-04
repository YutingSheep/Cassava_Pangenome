// i18n.ts - 国际化支持文件
import { create } from 'zustand';

// 定义语言类型
export type Language = 'zh' | 'en';

// 定义国际化存储接口
interface I18nStore {
  language: Language;
  setLanguage: (language: Language) => void;
}

// 创建国际化状态存储
export const useI18nStore = create<I18nStore>((set) => ({
  language: 'zh', // 默认语言为中文
  setLanguage: (language) => set({ language }),
}));

// 翻译数据
export const translations = {
  zh: {
    // 导航
    nav: {
      home: '首页',
      github: 'GitHub',
      documentation: '文档',
    },
    // 标题和简介
    hero: {
      title: '木薯泛基因组分析工具集',
      subtitle: '用于木薯泛基因组研究的综合分析脚本和工具',
      description: '本仓库收集了木薯泛基因组项目中使用的各类分析脚本和工具。木薯作为全球重要的粮食和工业原料作物，其基因组研究对于品种改良和产量提升具有重要意义。本项目旨在通过泛基因组分析，深入了解木薯种质资源的遗传多样性，为木薯育种和基因功能研究提供支持。',
      getStarted: '开始使用',
      viewOnGithub: '在GitHub上查看',
    },
    // 功能模块
    modules: {
      title: '主要功能模块',
      // 第一部分
      part1: {
        title: '基因组分析基础工具',
        module1: {
          title: '基因组组装',
          description: '高通量测序数据预处理、从头组装流程、混合组装策略和组装结果整合。'
        },
        module2: {
          title: '基因组质量评估',
          description: 'BUSCO完整性分析、基因组统计指标计算、组装连续性评估和测序深度分析。'
        },
        module3: {
          title: '基因组注释',
          description: '基因结构预测、功能注释流程、非编码RNA识别和重复序列注释。'
        },
        module4: {
          title: '泛基因组构建',
          description: '多样本比对策略、核心/可变基因组划分、基因簇分析和泛基因组可视化。'
        },
        module5: {
          title: '变异检测',
          description: 'SNP/InDel鉴定、结构变异分析、拷贝数变异检测和变异注释与过滤。'
        }
      },
      // 第二部分
      part2: {
        title: '群体遗传分析工具',
        module1: {
          title: '亲缘关系构建',
          description: '遗传距离计算、系统发育树构建、主成分分析和群体结构分析。'
        },
        module2: {
          title: '单倍型共享分析',
          description: '单倍型区块识别、连锁不平衡分析、单倍型网络构建和共享区段可视化。'
        },
        module3: {
          title: '杂合度计算',
          description: '个体杂合率估算、群体杂合度分布、杂合区域识别和纯合区段分析。'
        }
      },
      // 第三部分
      part3: {
        title: '进化分析工具',
        module1: {
          title: '全基因组复制(WGD)分析',
          description: '同源基因识别、复制事件定年、共线性区块分析和复制模式推断。'
        },
        module2: {
          title: '基因保留计算',
          description: '基因家族分析、保留/丢失模式识别、功能偏好性分析和进化速率计算。'
        },
        module3: {
          title: '选择压力分析',
          description: '正选择/负选择检测、适应性进化分析、基因功能富集和选择信号可视化。'
        }
      },
      // 第四部分
      part4: {
        title: '功能基因组学工具',
        module1: {
          title: '共线性分析',
          description: '种内/种间共线性比较、基因组重排检测、共线性区块可视化和进化断点分析。'
        },
        module2: {
          title: 'GWAS分析',
          description: '表型数据处理、关联分析流程、显著位点鉴定和候选基因筛选。'
        },
        module3: {
          title: '遗传多样性分析',
          description: '多样性参数计算、群体分化分析、选择清除检测和基因流分析。'
        }
      }
    },
    // 使用指南
    usage: {
      title: '使用指南',
      description: '每个模块下包含详细的使用说明文档，包括依赖安装、参数设置和结果解释。用户可根据研究需求选择相应的分析流程。',
      readDocs: '阅读文档'
    },
    // 贡献指南
    contribution: {
      title: '贡献指南',
      description: '欢迎研究人员提交改进建议或新的分析脚本。请遵循项目的代码规范和文档标准。',
      contribute: '如何贡献'
    },
    // 引用
    citation: {
      title: '引用',
      description: '如果您在研究中使用了本仓库的工具，请引用我们的相关论文。'
    },
    // 许可证
    license: {
      title: '许可证',
      description: '本项目采用 MIT 许可证。'
    },
    // 页脚
    footer: {
      copyright: '© 2025 木薯泛基因组项目。保留所有权利。'
    }
  },
  en: {
    // Navigation
    nav: {
      home: 'Home',
      github: 'GitHub',
      documentation: 'Documentation',
    },
    // Title and introduction
    hero: {
      title: 'Cassava Pangenome Analysis Toolkit',
      subtitle: 'Comprehensive analysis scripts and tools for cassava pangenome research',
      description: 'This repository contains a collection of scripts and tools used in the Cassava Pangenome Project. As a globally important food and industrial crop, cassava genome research is crucial for variety improvement and yield enhancement. This project aims to comprehensively understand the genetic diversity of cassava germplasm resources through pangenome analysis, providing support for cassava breeding and gene function research.',
      getStarted: 'Get Started',
      viewOnGithub: 'View on GitHub',
    },
    // Functional modules
    modules: {
      title: 'Main Functional Modules',
      // Part I
      part1: {
        title: 'Genome Analysis Fundamentals',
        module1: {
          title: 'Genome Assembly',
          description: 'High-throughput sequencing data preprocessing, de novo assembly pipeline, hybrid assembly strategies, and assembly integration.'
        },
        module2: {
          title: 'Genome Quality Assessment',
          description: 'BUSCO completeness analysis, genome statistics calculation, assembly continuity evaluation, and sequencing depth analysis.'
        },
        module3: {
          title: 'Genome Annotation',
          description: 'Gene structure prediction, functional annotation pipeline, non-coding RNA identification, and repetitive sequence annotation.'
        },
        module4: {
          title: 'Pangenome Construction',
          description: 'Multi-sample alignment strategies, core/variable genome partitioning, gene cluster analysis, and pangenome visualization.'
        },
        module5: {
          title: 'Variant Detection',
          description: 'SNP/InDel identification, structural variation analysis, copy number variation detection, and variant annotation and filtering.'
        }
      },
      // Part II
      part2: {
        title: 'Population Genetics Tools',
        module1: {
          title: 'Kinship Construction',
          description: 'Genetic distance calculation, phylogenetic tree construction, principal component analysis, and population structure analysis.'
        },
        module2: {
          title: 'Haplotype Sharing Analysis',
          description: 'Haplotype block identification, linkage disequilibrium analysis, haplotype network construction, and shared segment visualization.'
        },
        module3: {
          title: 'Heterozygosity Calculation',
          description: 'Individual heterozygosity estimation, population heterozygosity distribution, heterozygous region identification, and runs of homozygosity analysis.'
        }
      },
      // Part III
      part3: {
        title: 'Evolutionary Analysis Tools',
        module1: {
          title: 'Whole Genome Duplication (WGD) Analysis',
          description: 'Homologous gene identification, duplication event dating, syntenic block analysis, and duplication pattern inference.'
        },
        module2: {
          title: 'Gene Retention Calculation',
          description: 'Gene family analysis, retention/loss pattern identification, functional preference analysis, and evolutionary rate calculation.'
        },
        module3: {
          title: 'Selection Analysis',
          description: 'Positive/negative selection detection, adaptive evolution analysis, gene function enrichment, and selection signal visualization.'
        }
      },
      // Part IV
      part4: {
        title: 'Functional Genomics Tools',
        module1: {
          title: 'Synteny Analysis',
          description: 'Intra/inter-species synteny comparison, genome rearrangement detection, syntenic block visualization, and evolutionary breakpoint analysis.'
        },
        module2: {
          title: 'GWAS Analysis',
          description: 'Phenotype data processing, association analysis pipeline, significant locus identification, and candidate gene screening.'
        },
        module3: {
          title: 'Genetic Diversity Analysis',
          description: 'Diversity parameter calculation, population differentiation analysis, selective sweep detection, and gene flow analysis.'
        }
      }
    },
    // Usage guide
    usage: {
      title: 'Usage Guide',
      description: 'Each module includes detailed documentation covering dependency installation, parameter settings, and result interpretation. Users can select appropriate analysis workflows based on their research needs.',
      readDocs: 'Read Documentation'
    },
    // Contribution guidelines
    contribution: {
      title: 'Contribution Guidelines',
      description: 'Researchers are welcome to submit improvement suggestions or new analysis scripts. Please follow the project\'s code standards and documentation guidelines.',
      contribute: 'How to Contribute'
    },
    // Citation
    citation: {
      title: 'Citation',
      description: 'If you use tools from this repository in your research, please cite our relevant papers.'
    },
    // License
    license: {
      title: 'License',
      description: 'This project is licensed under the MIT License.'
    },
    // Footer
    footer: {
      copyright: '© 2025 Cassava Pangenome Project. All rights reserved.'
    }
  }
};
